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
	if (ops == R_NilValue) {
		for (int i = 0; i < 256; i++)
			ops_lkup_table[i] = 1;
		return;
	}
	for (int i = 0; i < 256; i++)
		ops_lkup_table[i] = 0;
	int ops_len = LENGTH(ops);
	for (int i = 0; i < ops_len; i++) {
		SEXP ops_elt = STRING_ELT(ops, i);
		if (ops_elt == NA_STRING || LENGTH(ops_elt) == 0)
			error("'ops' contains NAs and/or empty strings");
		char OP = CHAR(ops_elt)[0];
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
     cigars: character vector containing the extended CIGAR string for each
             read;
     ans_type: a single integer specifying the type of answer to return:
       0: 'ans' is a string describing the first validity failure or NULL;
       1: 'ans' is logical vector with TRUE values for valid elements
          in 'cigars'. */
SEXP C_validate_cigars(SEXP cigars, SEXP ans_type)
{
	SEXP ans;
	int ncigars, ans_type0, i;
	const char *cigar_string, *errmsg;
	char string_buf[200];

	ncigars = LENGTH(cigars);
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0 == 1)
		PROTECT(ans = NEW_LOGICAL(ncigars));
	else
		ans = R_NilValue;
	for (i = 0; i < ncigars; i++) {
		SEXP cigars_elt = STRING_ELT(cigars, i);
		if (cigars_elt == NA_STRING) {
			if (ans_type0 == 1)
				LOGICAL(ans)[i] = 1;
			continue;
		}
		cigar_string = CHAR(cigars_elt);
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
 *   cigars: character vector containing the extended CIGAR strings to
 *           explode.
 *   ops:    NULL or a character vector containing the CIGAR operations to
 *           actually consider. If NULL, then all CIGAR operations are
 *           considered.
 * Both functions return a list of the same length as 'cigars' where each
 * list element is a character vector (for C_explode_cigar_ops()) or an integer
 * vector (for C_explode_cigar_oplens()). The 2 lists have the same shape,
 * that is, same length() and same elementNROWS(). The i-th character vector
 * in the list returned by C_explode_cigar_ops() contains one single-letter
 * string per CIGAR operation in 'cigars[i]'. The i-th integer vector in the
 * list returned by C_explode_cigar_oplens() contains the corresponding
 * CIGAR operation lengths. Zero-length operations or operations not listed
 * in 'ops' are ignored.
 */
SEXP C_explode_cigar_ops(SEXP cigars, SEXP ops)
{
	SEXP ans, ans_elt, ans_elt_elt;
	int ncigars, ans_elt_len, i, j;
	CharAE *OPbuf;
	const char *cigar_string, *errmsg;

	ncigars = LENGTH(cigars);
	_init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(ncigars));
	OPbuf = new_CharAE(0);
	for (i = 0; i < ncigars; i++) {
		SEXP cigars_elt = STRING_ELT(cigars, i);
		if (cigars_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigars[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigars_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigars[%d]' is \"*\"", i + 1);
		}
		CharAE_set_nelt(OPbuf, 0);
		errmsg = split_cigar_string(cigar_string, OPbuf, NULL);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
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

SEXP C_explode_cigar_oplens(SEXP cigars, SEXP ops)
{
	SEXP ans, ans_elt;
	int ncigars, i;
	IntAE *OPLbuf;
	const char *cigar_string, *errmsg;

	ncigars = LENGTH(cigars);
	_init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(ncigars));
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < ncigars; i++) {
		SEXP cigars_elt = STRING_ELT(cigars, i);
		if (cigars_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigars[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigars_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigars[%d]' is \"*\"", i + 1);
		}
		IntAE_set_nelt(OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, NULL, OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
		}
		PROTECT(ans_elt = new_INTEGER_from_IntAE(OPLbuf));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

