#include "tabulate_cigar_ops.h"

#include "inspect_cigars.h"

#include <string.h>  /* for memset() */


static const char *cigar_string_op_table(SEXP cigar_string, int weighted,
		const char *allOPs, int *table_row, int table_nrow)
{
	static char errmsg_buf[200];

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	const char *cig0 = CHAR(cigar_string);
	char OP /* Operation */;
	int n, offset = 0, OPL /* Operation Length */;
	while ((n = _next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		const char *tmp = strchr(allOPs, (int) OP);
		if (tmp == NULL) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		int *dest_p = table_row + (tmp - allOPs) * table_nrow;
		*dest_p += weighted ? OPL : 1;
		offset += n;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigars: character vector containing the extended CIGAR string for each
 *           read;
 * Return an integer matrix with the number of rows equal to the length of
 * 'cigars' and 9 columns, one for each extended CIGAR operation containing
 * a frequency count for the operations for each element of 'cigars'.
 */
SEXP C_tabulate_cigar_ops(SEXP cigars, SEXP oplens_as_weights)
{
	static const char *allOPs = "MIDNSHP=X";

	int cigar_len = LENGTH(cigars);
	int weighted = LOGICAL(oplens_as_weights)[0];
	int allOPs_len = strlen(allOPs);
	SEXP ans = PROTECT(allocMatrix(INTSXP, cigar_len, allOPs_len));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	int *ans_p = INTEGER(ans);
	for (int i = 0; i < cigar_len; i++, ans_p++) {
		SEXP cigar_string = STRING_ELT(cigars, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		const char *errmsg = cigar_string_op_table(cigar_string,
							   weighted, allOPs,
							   ans_p, cigar_len);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
		}
	}

	SEXP ans_colnames = PROTECT(NEW_CHARACTER(allOPs_len));
	for (int j = 0; j < allOPs_len; j++) {
		SEXP OP = PROTECT(mkCharLen(allOPs + j, 1));
		SET_STRING_ELT(ans_colnames, j, OP);
		UNPROTECT(1);
	}
	SEXP ans_dimnames = PROTECT(NEW_LIST(2));
	SET_ELEMENT(ans_dimnames, 0, R_NilValue);
	SET_ELEMENT(ans_dimnames, 1, ans_colnames);
	SET_DIMNAMES(ans, ans_dimnames);
	UNPROTECT(3);
	return ans;
}

