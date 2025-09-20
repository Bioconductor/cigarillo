#include "trim_cigars.h"

#include "inspect_cigars.h"


static char errmsg_buf[200];


/****************************************************************************
 * C_trim_cigars_along_ref()
 */

static const char *Ltrim_along_ref(SEXP cigar_string,
				   int *Lnpos, int *Loffset, int *rshift)
{
	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	const char *cig0 = CHAR(cigar_string);
	*rshift = 0;
	int n, offset = 0, OPL /* Operation Length */;
	char OP /* Operation */;
	while ((n = _next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Lnpos < OPL) {
				*Loffset = offset;
				*rshift += *Lnpos;
				return NULL;
			}
			*Lnpos -= OPL;
			*rshift += OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			if (*Lnpos < OPL)
				*Lnpos = 0;
			else
				*Lnpos -= OPL;
			*rshift += OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		offset += n;
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after trimming");
	return errmsg_buf;
}

static const char *Rtrim_along_ref(SEXP cigar_string,
				   int *Rnpos, int *Roffset)
{
	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	const char *cig0 = CHAR(cigar_string);
	int n, offset = LENGTH(cigar_string), OPL /* Operation Length */;
	char OP /* Operation */;
	while ((n = _prev_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		offset -= n;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Rnpos < OPL) {
				*Roffset = offset;
				return NULL;
			}
			*Rnpos -= OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			if (*Rnpos < OPL)
				*Rnpos = 0;
			else
				*Rnpos -= OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after trimming");
	return errmsg_buf;
}

#define	CIGAR_BUF_LENGTH 250000

/* FIXME: 'cigar_buf' is under the risk of a buffer overflow! */
static const char *trim_along_ref(SEXP cigar_string,
		int Lnpos, int Rnpos, char *cigar_buf, int *rshift)
{
	//Rprintf("trim_along_ref():\n");
	int Loffset;
	const char *errmsg =
		Ltrim_along_ref(cigar_string, &Lnpos, &Loffset, rshift);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Lnpos=%d Loffset=%d *rshift=%d\n",
	//	Lnpos, Loffset, *rshift);
	int Roffset;
	errmsg = Rtrim_along_ref(cigar_string, &Rnpos, &Roffset);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Rnpos=%d Roffset=%d\n", Rnpos, Roffset);
	if (Roffset < Loffset) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "CIGAR is empty after trimming");
		return errmsg_buf;
	}
	int buf_offset = 0, n;
	const char *cig0 = CHAR(cigar_string);
	for (int offset = Loffset; offset <= Roffset; offset += n) {
		char OP /* Operation */;
		int OPL /* Operation Length */;
		n = _next_cigar_OP(cig0, offset, &OP, &OPL);
		if (offset == Loffset)
			OPL -= Lnpos;
		if (offset == Roffset)
			OPL -= Rnpos;
		if (OPL <= 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "CIGAR is empty after trimming");
			return errmsg_buf;
		}
		size_t size = CIGAR_BUF_LENGTH - buf_offset;
		int ret = snprintf(cigar_buf + buf_offset, size,
				   "%d%c", OPL, OP);
		if (ret >= size) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'cigar_buf' overflow");
			return errmsg_buf;
		}
		buf_offset += ret;
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP C_trim_cigars_along_ref(SEXP cigars, SEXP Lnpos, SEXP Rnpos)
{
	static char cigar_buf[CIGAR_BUF_LENGTH];

	int ncigars = LENGTH(cigars);
	SEXP trimmed_cigars = PROTECT(NEW_CHARACTER(ncigars));
	SEXP ans_rshift = PROTECT(NEW_INTEGER(ncigars));
	for (int i = 0; i < ncigars; i++) {
		SEXP cigar_string = STRING_ELT(cigars, i);
		if (cigar_string == NA_STRING) {
			SET_STRING_ELT(trimmed_cigars, i, NA_STRING);
			INTEGER(ans_rshift)[i] = NA_INTEGER;
			continue;
		}
		const char *errmsg = trim_along_ref(cigar_string,
					INTEGER(Lnpos)[i],
					INTEGER(Rnpos)[i],
					cigar_buf, INTEGER(ans_rshift) + i);
		if (errmsg != NULL) {
			UNPROTECT(2);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
		}
		SEXP trimmed_string = PROTECT(mkChar(cigar_buf));
		SET_STRING_ELT(trimmed_cigars, i, trimmed_string);
		UNPROTECT(1);
	}

	SEXP ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, trimmed_cigars);
	SET_VECTOR_ELT(ans, 1, ans_rshift);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * C_trim_cigars_along_query()
 */

static const char *Ltrim_along_query(SEXP cigar_string,
				     int *Lnpos, int *Loffset, int *rshift)
{
	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	const char *cig0 = CHAR(cigar_string);
	*rshift = 0;
	int n, offset = 0, OPL /* Operation Length */;
	char OP /* Operation */;
	while ((n = _next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Lnpos < OPL) {
				*Loffset = offset;
				*rshift += *Lnpos;
				return NULL;
			}
			*Lnpos -= OPL;
			*rshift += OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
			if (*Lnpos < OPL) {
				*Loffset = offset;
				return NULL;
			}
			*Lnpos -= OPL;
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			*rshift += OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		offset += n;
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after trimming");
	return errmsg_buf;
}

static const char *Rtrim_along_query(SEXP cigar_string,
				     int *Rnpos, int *Roffset)
{
	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	const char *cig0 = CHAR(cigar_string);
	int n, offset = LENGTH(cigar_string), OPL /* Operation Length */;
	char OP /* Operation */;
	while ((n = _prev_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		offset -= n;
		switch (OP) {
		/* M, =, X, I, S, H */
		    case 'M': case '=': case 'X': case 'I': case 'S': case 'H':
			if (*Rnpos < OPL) {
				*Roffset = offset;
				return NULL;
			}
			*Rnpos -= OPL;
		    break;
		/* Deletion (or skipped region) from the reference,
		   or silent deletion from the padded reference */
		    case 'D': case 'N': case 'P':
		    break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after trimming");
	return errmsg_buf;
}

/* FIXME: 'cigar_buf' is under the risk of a buffer overflow! */
static const char *trim_along_query(SEXP cigar_string,
		int Lnpos, int Rnpos, char *cigar_buf, int *rshift)
{
	//Rprintf("trim_along_query():\n");
	int Loffset;
	const char *errmsg =
		Ltrim_along_query(cigar_string, &Lnpos, &Loffset, rshift);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Lnpos=%d Loffset=%d *rshift=%d\n",
	//	Lnpos, Loffset, *rshift);
	int Roffset;
	errmsg = Rtrim_along_query(cigar_string, &Rnpos, &Roffset);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Rnpos=%d Roffset=%d\n", Rnpos, Roffset);
	if (Roffset < Loffset) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "CIGAR is empty after trimming");
		return errmsg_buf;
	}
	int buf_offset = 0, n;
	const char *cig0 = CHAR(cigar_string);
	for (int offset = Loffset; offset <= Roffset; offset += n) {
		char OP /* Operation */;
		int OPL /* Operation Length */;
		n = _next_cigar_OP(cig0, offset, &OP, &OPL);
		if (offset == Loffset)
			OPL -= Lnpos;
		if (offset == Roffset)
			OPL -= Rnpos;
		if (OPL <= 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "CIGAR is empty after trimming");
			return errmsg_buf;
		}
		size_t size = CIGAR_BUF_LENGTH - buf_offset;
		int ret = snprintf(cigar_buf + buf_offset, size,
				   "%d%c", OPL, OP);
		if (ret >= size) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'cigar_buf' overflow");
			return errmsg_buf;
		}
		buf_offset += ret;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
   Return a list of 2 elements:
     1. The vector of trimmed CIGARs.
     2. The 'rshift' vector i.e. the integer vector of the same length
        as 'cigars' that would need to be added to the 'lmmpos' field
        of a SAM/BAM file as a consequence of this trimming. */
SEXP C_trim_cigars_along_query(SEXP cigars, SEXP Lnpos, SEXP Rnpos)
{
	static char cigar_buf[CIGAR_BUF_LENGTH];

	int ncigars = LENGTH(cigars);
	SEXP trimmed_cigars = PROTECT(NEW_CHARACTER(ncigars));
	SEXP ans_rshift = PROTECT(NEW_INTEGER(ncigars));
	for (int i = 0; i < ncigars; i++) {
		SEXP cigar_string = STRING_ELT(cigars, i);
		if (cigar_string == NA_STRING) {
			SET_STRING_ELT(trimmed_cigars, i, NA_STRING);
			INTEGER(ans_rshift)[i] = NA_INTEGER;
			continue;
		}
		const char *errmsg = trim_along_query(cigar_string,
					INTEGER(Lnpos)[i],
					INTEGER(Rnpos)[i],
					cigar_buf, INTEGER(ans_rshift) + i);
		if (errmsg != NULL) {
			UNPROTECT(2);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
		}
		SEXP trimmed_string = PROTECT(mkChar(cigar_buf));
		SET_STRING_ELT(trimmed_cigars, i, trimmed_string);
		UNPROTECT(1);
	}

	SEXP ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, trimmed_cigars);
	SET_VECTOR_ELT(ans, 1, ans_rshift);
	UNPROTECT(3);
	return ans;
}

