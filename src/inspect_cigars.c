#include "inspect_cigars.h"

#include "S4Vectors_interface.h"

#include <ctype.h> /* for isdigit() */


static char errmsg_buf[200];


const char *_get_cigar_parsing_error()
{
	return errmsg_buf;
}

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. cigar_string[offset] is '\0'), or -1 in case of a parse error
   (in which case _get_cigar_parsing_error() can be used to get a pointer to
   the error message).
   Zero-length operations are ignored. */
int _next_cigar_OP(const char *cigar_string, int offset, char *OP, int *OPL)
{
	char c;
	int offset0, opl;

	if (!cigar_string[offset])
		return 0;
	offset0 = offset;
	do {
		/* Extract *OPL */
		opl = 0;
		while (isdigit(c = cigar_string[offset])) {
			offset++;
			opl *= 10;
			opl += c - '0';
		}
		/* Extract *OP */
		if (!(*OP = cigar_string[offset])) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unexpected CIGAR end after char %d",
				 offset);
			return -1;
		}
		offset++;
	} while (opl == 0);
	*OPL = opl;
	return offset - offset0;
}

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. offset is 0), or -1 in case of a parse error.
   Zero-length operations are ignored. */
int _prev_cigar_OP(const char *cigar_string, int offset, char *OP, int *OPL)
{
	char c;
	int offset0, opl, powof10;

	if (offset == 0)
		return 0;
	offset0 = offset;
	do {
		/* Extract *OP */
		offset--;
		*OP = cigar_string[offset];
		/* Extract *OPL */
		if (offset == 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "no CIGAR operation length before char %d",
				 offset + 1);
			return -1;
		}
		offset--;
		opl = 0;
		powof10 = 1;
		while (offset >= 0 && isdigit(c = cigar_string[offset])) {
			opl += (c - '0') * powof10;
			powof10 *= 10;
			offset--;
		}
		offset++;
	} while (opl == 0);
	*OPL = opl;
	return offset0 - offset;
}


/****************************************************************************
 * _is_in_ops()
 */

static int ops_lkup_table[256];

void _init_ops_lkup_table(SEXP ops)
{
	int ops_len, i;
	SEXP ops_elt;
	char OP;

	if (ops == R_NilValue) {
		for (i = 0; i < 256; i++)
			ops_lkup_table[i] = 1;
		return;
	}
	for (i = 0; i < 256; i++)
		ops_lkup_table[i] = 0;
	ops_len = LENGTH(ops);
	for (i = 0; i < ops_len; i++) {
		ops_elt = STRING_ELT(ops, i);
		if (ops_elt == NA_STRING || LENGTH(ops_elt) == 0)
			error("'ops' contains NAs and/or empty strings");
		OP = CHAR(ops_elt)[0];
		ops_lkup_table[(unsigned char) OP] = 1;
	}
	return;
}

/* A fast way to determine whether a CIGAR operation is in the 'ops' vector
   that was preprocessed with '_init_ops_lkup_table(ops)'. So make sure to
   initialize 'ops_lkup_table' with '_init_ops_lkup_table(ops)' before
   calling '_is_in_ops(OP)'. */
int _is_in_ops(char OP)
{
	return ops_lkup_table[(unsigned char) OP];
}


/****************************************************************************
 * C_validate_cigars()
 */

static const char *parse_cigar(const char *cigar_string)
{
	int cigar_offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	cigar_offset = 0;
	while ((n = _next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
   Args:
     cigar: character vector containing the extended CIGAR string for each
            read;
     ans_type: a single integer specifying the type of answer to return:
       0: 'ans' is a string describing the first validity failure or NULL;
       1: 'ans' is logical vector with TRUE values for valid elements
          in 'cigar'. */
SEXP C_validate_cigars(SEXP cigar, SEXP ans_type)
{
	SEXP ans, cigar_elt;
	int cigar_len, ans_type0, i;
	const char *cigar_string, *errmsg;
	char string_buf[200];

	cigar_len = LENGTH(cigar);
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0 == 1)
		PROTECT(ans = NEW_LOGICAL(cigar_len));
	else
		ans = R_NilValue;
	for (i = 0; i < cigar_len; i++) {
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			if (ans_type0 == 1)
				LOGICAL(ans)[i] = 1;
			continue;
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			if (ans_type0 == 1)
				LOGICAL(ans)[i] = 1;
			continue;
		}
		errmsg = parse_cigar(cigar_string);
		if (ans_type0 == 1) {
			LOGICAL(ans)[i] = errmsg == NULL;
			continue;
		}
		if (errmsg != NULL) {
			snprintf(string_buf, sizeof(string_buf),
				 "element %d is invalid (%s)", i + 1, errmsg);
			return mkString(string_buf);
		}
	}
	if (ans_type0 == 1)
		UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_explode_cigar_ops() and C_explode_cigar_oplens()
 */

/* Make sure _init_ops_lkup_table() is called before split_cigar_string(). */
static const char *split_cigar_string(const char *cigar_string,
		CharAE *OPbuf, IntAE *OPLbuf)
{
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	offset = 0;
	while ((n = _next_cigar_OP(cigar_string, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		if (_is_in_ops(OP)) {
			if (OPbuf != NULL)
				CharAE_insert_at(OPbuf,
					CharAE_get_nelt(OPbuf), OP);
			if (OPLbuf != NULL)
				IntAE_insert_at(OPLbuf,
					IntAE_get_nelt(OPLbuf), OPL);
		}
		offset += n;
	}
	return NULL;
}

/* --- .Call ENTRY POINTS ---
 *   - C_explode_cigar_ops()
 *   - C_explode_cigar_oplens()
 * Args:
 *   cigar: character vector containing the extended CIGAR strings to
 *          explode.
 *   ops:   NULL or a character vector containing the CIGAR operations to
 *          actually consider. If NULL, then all CIGAR operations are
 *          considered.
 * Both functions return a list of the same length as 'cigar' where each
 * list element is a character vector (for C_explode_cigar_ops()) or an integer
 * vector (for C_explode_cigar_oplens()). The 2 lists have the same shape,
 * that is, same length() and same elementNROWS(). The i-th character vector
 * in the list returned by C_explode_cigar_ops() contains one single-letter
 * string per CIGAR operation in 'cigar[i]'. The i-th integer vector in the
 * list returned by C_explode_cigar_oplens() contains the corresponding
 * CIGAR operation lengths. Zero-length operations or operations not listed
 * in 'ops' are ignored.
 */
SEXP C_explode_cigar_ops(SEXP cigar, SEXP ops)
{
	SEXP ans, cigar_elt, ans_elt, ans_elt_elt;
	int cigar_len, ans_elt_len, i, j;
	CharAE *OPbuf;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	_init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(cigar_len));
	OPbuf = new_CharAE(0);
	for (i = 0; i < cigar_len; i++) {
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		CharAE_set_nelt(OPbuf, 0);
		errmsg = split_cigar_string(cigar_string, OPbuf, NULL);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
		ans_elt_len = CharAE_get_nelt(OPbuf);
		PROTECT(ans_elt = NEW_CHARACTER(ans_elt_len));
		for (j = 0; j < ans_elt_len; j++) {
			PROTECT(ans_elt_elt = mkCharLen(OPbuf->elts + j, 1));
			SET_STRING_ELT(ans_elt, j, ans_elt_elt);
			UNPROTECT(1);
		}
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

SEXP C_explode_cigar_oplens(SEXP cigar, SEXP ops)
{
	SEXP ans, cigar_elt, ans_elt;
	int cigar_len, i;
	IntAE *OPLbuf;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	_init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(cigar_len));
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_len; i++) {
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		IntAE_set_nelt(OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, NULL, OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
		PROTECT(ans_elt = new_INTEGER_from_IntAE(OPLbuf));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_tabulate_cigar_ops()
 */

static const char *cigar_string_op_table(SEXP cigar_string, const char *allOPs,
		int *table_row, int table_nrow)
{
	const char *cig0, *tmp;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	offset = 0;
	while ((n = _next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		tmp = strchr(allOPs, (int) OP);
		if (tmp == NULL) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		*(table_row + (tmp - allOPs) * table_nrow) += OPL;
		offset += n;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 * Return an integer matrix with the number of rows equal to the length of
 * 'cigar' and 9 columns, one for each extended CIGAR operation containing
 * a frequency count for the operations for each element of 'cigar'.
 */
SEXP C_tabulate_cigar_ops(SEXP cigar)
{
	SEXP cigar_string, ans, ans_dimnames, ans_colnames;
	int cigar_len, allOPs_len, i, j, *ans_row;
	const char *allOPs = "MIDNSHP=X", *errmsg;
	char OPstrbuf[2];

	cigar_len = LENGTH(cigar);
	allOPs_len = strlen(allOPs);
	PROTECT(ans = allocMatrix(INTSXP, cigar_len, allOPs_len));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	ans_row = INTEGER(ans);
	for (i = 0, ans_row = INTEGER(ans); i < cigar_len; i++, ans_row++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_op_table(cigar_string, allOPs,
				ans_row, cigar_len);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
	}

	PROTECT(ans_colnames = NEW_CHARACTER(allOPs_len));
	OPstrbuf[1] = '\0';
	for (j = 0; j < allOPs_len; j++) {
		OPstrbuf[0] = allOPs[j];
		SET_STRING_ELT(ans_colnames, j, mkChar(OPstrbuf));
	}
	PROTECT(ans_dimnames = NEW_LIST(2));
	SET_ELEMENT(ans_dimnames, 0, R_NilValue);
	SET_ELEMENT(ans_dimnames, 1, ans_colnames);
	SET_DIMNAMES(ans, ans_dimnames);
	UNPROTECT(3);
	return ans;
}

