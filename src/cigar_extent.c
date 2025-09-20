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
   cigars, space, flags: see C_cigars_as_ranges() in src/cigars_as_ranges.c
   Returns an integer vector of the same length as 'cigars' containing the
   extents of the alignments as inferred from the cigars information. */
SEXP C_cigar_extent(SEXP cigars, SEXP space, SEXP flags)
{
	SEXP ans;
	int space0, i, *ans_elt;
	const int *flags_elt;
	const char *cigar_string, *errmsg;

	int ncigars = LENGTH(cigars);
	if (flags != R_NilValue)
		flags_elt = INTEGER(flags);
	space0 = INTEGER(space)[0];
	PROTECT(ans = NEW_INTEGER(ncigars));
	for (i = 0, ans_elt = INTEGER(ans); i < ncigars; i++, ans_elt++) {
		if (flags != R_NilValue) {
			if (*flags_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flags' contains NAs");
			}
			if (*flags_elt & 0x004) {
				*ans_elt = NA_INTEGER;
				goto for_tail;
			}
		}
		SEXP cigars_elt = STRING_ELT(cigars, i);
		if (cigars_elt == NA_STRING) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		cigar_string = CHAR(cigars_elt);
		if (strcmp(cigar_string, "*") == 0) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		errmsg = compute_cigar_extent(cigar_string, space0, ans_elt);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigars[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flags != R_NilValue)
			flags_elt++;
	}
	UNPROTECT(1);
	return ans;
}

