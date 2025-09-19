#include "cigar_extent.h"

#include "cigar_ops_visibility.h"
#include "inspect_cigars.h"


static const char *compute_cigar_extent(const char *cigar_string, int space,
					int *extent)
{
	int x, cigar_offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	x = cigar_offset = 0;
	while ((n = _next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		if (_op_is_visible(OP, space))
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

