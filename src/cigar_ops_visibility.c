#include "cigar_ops_visibility.h"


/* The 8 "projection spaces" below are also defined at the top of the
   R/project_sequences.R file. */
#define REFERENCE                       1
#define REFERENCE_N_REGIONS_REMOVED     2
#define QUERY                           3
#define QUERY_BEFORE_HARD_CLIPPING      4
#define QUERY_AFTER_SOFT_CLIPPING       5
#define PAIRWISE                        6
#define PAIRWISE_N_REGIONS_REMOVED      7
#define PAIRWISE_DENSE                  8


/****************************************************************************
 * C_cigar_ops_visibility()
 */

int _op_is_visible(char OP, int space)
{
	if (OP == 'M')
		return 1;
	switch (space) {
	    case QUERY_BEFORE_HARD_CLIPPING:
		if (OP == 'H')
			return 1;
		/* fall through */
	    case QUERY:
		if (OP == 'S')
			return 1;
		/* fall through */
	    case QUERY_AFTER_SOFT_CLIPPING:
		if (OP == 'I')
			return 1;
		break;
	    case PAIRWISE:
		if (OP == 'I')
			return 1;
		/* fall through */
	    case REFERENCE:
		if (OP == 'D' || OP == 'N')
			return 1;
		break;
	    case PAIRWISE_N_REGIONS_REMOVED:
		if (OP == 'I')
			return 1;
		/* fall through */
	    case REFERENCE_N_REGIONS_REMOVED:
		if (OP == 'D')
			return 1;
	}
	if (OP == '=' || OP == 'X')
		return 1;
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_cigar_ops_visibility(SEXP ops)
{
	int ops_len = LENGTH(ops);
	SEXP ans = PROTECT(allocMatrix(INTSXP, 8, ops_len));
	int *ans_p = INTEGER(ans);
	for (int j = 0; j < ops_len; j++) {
		SEXP ops_elt = STRING_ELT(ops, j);
		if (ops_elt == NA_STRING || LENGTH(ops_elt) == 0) {
			UNPROTECT(1);
			error("'ops' contains NAs and/or empty strings");
		}
		char OP = CHAR(ops_elt)[0];
		for (int space = 1; space <= 8; space++)
			*(ans_p++) = _op_is_visible(OP, space);
	}
	UNPROTECT(1);
	return ans;
}

