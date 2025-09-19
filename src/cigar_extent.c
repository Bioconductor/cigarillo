#include "cigar_extent.h"

#include "inspect_cigars.h"


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

int _is_visible_in_space(char OP, int space)
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
		for (int i = 0; i < 8; i++)
			*(ans_p++) = _is_visible_in_space(OP, i + 1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_cigar_extent()
 */

static const char *compute_cigar_extent(const char *cigar_string, int space,
					int *extent)
{
	int x, cigar_offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	x = cigar_offset = 0;
	while ((n = _next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		if (_is_visible_in_space(OP, space))
			x += OPL;
		cigar_offset += n;
	}
	*extent = x;
	return NULL;
}

/* --- .Call ENTRY POINT ---
   Args:
   cigar, space, flag: see C_cigars_as_ranges() in src/cigars_as_ranges.c
   Returns an integer vector of the same length as 'cigar' containing the
   extents of the alignments as inferred from the cigar information. */
SEXP C_cigar_extent(SEXP cigar, SEXP space, SEXP flag)
{
	SEXP ans, cigar_elt;
	int cigar_len, space0, i, *ans_elt;
	const int *flag_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	if (flag != R_NilValue)
		flag_elt = INTEGER(flag);
	space0 = INTEGER(space)[0];
	PROTECT(ans = NEW_INTEGER(cigar_len));
	for (i = 0, ans_elt = INTEGER(ans); i < cigar_len; i++, ans_elt++) {
		if (flag != R_NilValue) {
			if (*flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (*flag_elt & 0x004) {
				*ans_elt = NA_INTEGER;
				goto for_tail;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		errmsg = compute_cigar_extent(cigar_string, space0, ans_elt);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flag != R_NilValue)
			flag_elt++;
	}
	UNPROTECT(1);
	return ans;
}

