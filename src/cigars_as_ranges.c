#include "cigars_as_ranges.h"

#include "IRanges_interface.h"
#include "S4Vectors_interface.h"

#include "inspect_cigars.h"
#include "cigar_extent.h"


static void drop_or_append_or_merge_range(int start, int width,
		int drop_empty_range, int merge_range, int nelt0,
		IntPairAE *range_buf, const char *OP, CharAEAE *OP_buf)
{
	int buf_nelt, buf_nelt_minus_1, prev_end_plus_1;
	CharAE *OP_buf_new_elt, *OP_buf_prev_elt;

	if (drop_empty_range && width == 0)  /* Drop. */
		return;
	buf_nelt = IntPairAE_get_nelt(range_buf);
	if (merge_range && buf_nelt > nelt0) {
		/* The incoming range should never overlap with the previous
		   incoming range i.e. 'start' should always be > the end of
		   the previous incoming range. */
		buf_nelt_minus_1 = buf_nelt - 1;
		prev_end_plus_1 = range_buf->a->elts[buf_nelt_minus_1] +
				  range_buf->b->elts[buf_nelt_minus_1];
		if (start == prev_end_plus_1) {
			/* Merge. */
			range_buf->b->elts[buf_nelt_minus_1] += width;
			if (OP_buf != NULL) {
				OP_buf_prev_elt =
					OP_buf->elts[buf_nelt_minus_1];
				CharAE_insert_at(OP_buf_prev_elt,
					CharAE_get_nelt(OP_buf_prev_elt), *OP);
			}
			return;
		}
	}
	/* Append. */
	IntPairAE_insert_at(range_buf, buf_nelt, start, width);
	if (OP_buf != NULL) {
		OP_buf_new_elt = new_CharAE(1);
		CharAE_insert_at(OP_buf_new_elt, 0, *OP);
		CharAEAE_insert_at(OP_buf, buf_nelt, OP_buf_new_elt);
	}
	return;
}

/* Make sure _init_ops_lkup_table() is called before parse_cigar_ranges(). */
static const char *parse_cigar_ranges(const char *cigar_string,
		int space, int lmmpos,
		int drop_empty_ranges, int reduce_ranges,
		IntPairAE *range_buf, CharAEAE *OP_buf)
{
	int buf_nelt0, cigar_offset, n, OPL /* Operation Length */,
	    start, width;
	char OP /* Operation */;

	buf_nelt0 = IntPairAE_get_nelt(range_buf);
	cigar_offset = 0;
	start = lmmpos;
	while ((n = _next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		width = _is_visible_in_space(OP, space) ? OPL : 0;
		if (_is_in_ops(OP))
			drop_or_append_or_merge_range(start, width,
						      drop_empty_ranges,
						      reduce_ranges, buf_nelt0,
						      range_buf, &OP, OP_buf);
		start += width;
		cigar_offset += n;
	}
	return NULL;
}


/****************************************************************************
 * C_cigars_as_ranges()
 */

static SEXP make_list_of_IRanges(const IntPairAEAE *range_buf, SEXP names)
{
	SEXP ans, ans_names;

	PROTECT(ans = new_list_of_IRanges_from_IntPairAEAE("IRanges",
							   range_buf));
	PROTECT(ans_names = duplicate(names));
	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}

static SEXP make_CompressedIRangesList(const IntPairAE *range_buf,
		const CharAEAE *OP_buf, SEXP breakpoints)
{
	SEXP ans, ans_unlistData, ans_unlistData_names, ans_partitioning;

	PROTECT(ans_unlistData =
			new_IRanges_from_IntPairAE("IRanges", range_buf));
	if (OP_buf != NULL) {
		PROTECT(ans_unlistData_names =
				new_CHARACTER_from_CharAEAE(OP_buf));
		set_IRanges_names(ans_unlistData, ans_unlistData_names);
		UNPROTECT(1);
	}
	PROTECT(ans_partitioning =
			new_PartitioningByEnd("PartitioningByEnd",
					      breakpoints, NULL));
	PROTECT(ans = new_CompressedList(
				"CompressedIRangesList",
				ans_unlistData, ans_partitioning));
	UNPROTECT(3);
	return ans;
}

/* --- .Call ENTRY POINT ---
   Args:
     cigar:  character vector containing extended CIGAR strings.
     flag:   NULL or an integer vector of the same length as 'cigar'
             containing the SAM flag for each read. Serves only as a way to
             indicate whether a read is mapped or not. According to the SAM
             Spec v1.4, flag bit 0x4 is the only reliable place to tell
             whether a segment (or read) is mapped (bit is 0) or not (bit is 1).
     space:  single integer indicating one of the 8 supported spaces (defined
             at the top of this file).
     lmmpos: integer vector of the same length as 'cigar' (or of length 1)
             containing the 1-based leftmost position/coordinate of the
             clipped read sequences.
     f:      NULL or a factor of length 'cigar'. If NULL, then the ranges are
             grouped by alignment and stored in a CompressedIRangesList object
             with 1 list element per element in 'cigar'. If a factor, then they
             are grouped by factor level and stored in an ordinary list of
             IRanges objects with 1 list element per level in 'f' and named
             with those levels.
     ops:    NULL or a character vector containing the CIGAR operations to
             translate to ranges. If NULL, then all CIGAR operations are
             translated.
     drop_empty_ranges: TRUE or FALSE.
     reduce_ranges: TRUE or FALSE.
     with_ops: TRUE or FALSE indicating whether the returned ranges should be
             named with their corresponding CIGAR operation.

   Returns either a CompressedIRangesList object of the same length as 'cigar'
   (if 'f' is NULL) or an ordinary list of IRanges objects with 1 list element
   per level in 'f' (if 'f' is a factor). This list is then turned into a
   SimpleIRangesList object in R. */
SEXP C_cigars_as_ranges(SEXP cigar, SEXP flag, SEXP space, SEXP lmmpos, SEXP f,
		SEXP ops, SEXP drop_empty_ranges, SEXP reduce_ranges,
		SEXP with_ops)
{
	SEXP ans, ans_breakpoints, f_levels, cigar_elt;
	int cigar_len, space0, pos_len, f_is_NULL, ans_len, *breakpoint,
	    drop_empty_ranges0, reduce_ranges0, with_ops0, i;
	IntPairAE *range_buf1;
	IntPairAEAE *range_buf2;
	CharAEAE *OP_buf;
	const int *flag_elt, *pos_elt, *f_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	if (flag != R_NilValue)
		flag_elt = INTEGER(flag);
	_init_ops_lkup_table(ops);
	space0 = INTEGER(space)[0];
	pos_len = LENGTH(lmmpos);
	pos_elt = INTEGER(lmmpos);
	f_is_NULL = f == R_NilValue;
	if (f_is_NULL) {
		ans_len = cigar_len;
		/* We will typically generate at least 'cigar_len' ranges. */
		range_buf1 = new_IntPairAE(ans_len, 0);
		PROTECT(ans_breakpoints = NEW_INTEGER(ans_len));
		breakpoint = INTEGER(ans_breakpoints);
	} else {
		f_levels = GET_LEVELS(f);
		ans_len = LENGTH(f_levels);
		range_buf2 = new_IntPairAEAE(ans_len, ans_len);
		f_elt = INTEGER(f);
	}
	drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	with_ops0 = LOGICAL(with_ops)[0];
	if (with_ops0 && f_is_NULL) {
		OP_buf = new_CharAEAE(cigar_len, 0);
	} else {
		OP_buf = NULL;
	}
	for (i = 0; i < cigar_len; i++) {
		if (flag != R_NilValue) {
			if (*flag_elt == NA_INTEGER) {
				if (f_is_NULL)
					UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (*flag_elt & 0x004) {
				/* The CIGAR of an unmapped read doesn't
				   produce any range i.e. it's treated as an
				   empty CIGAR. */
				goto for_tail;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		if (*pos_elt == NA_INTEGER || *pos_elt == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'lmmpos[%d]' is NA or 0", i + 1);
		}
		if (!f_is_NULL) {
			if (*f_elt == NA_INTEGER)
				error("'f[%d]' is NA", i + 1);
			range_buf1 = range_buf2->elts[*f_elt - 1];
		}
		errmsg = parse_cigar_ranges(cigar_string, space0, *pos_elt,
				drop_empty_ranges0, reduce_ranges0,
				range_buf1, OP_buf);
		if (errmsg != NULL) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flag != R_NilValue)
			flag_elt++;
		if (pos_len != 1)
			pos_elt++;
		if (f_is_NULL)
			*(breakpoint++) = IntPairAE_get_nelt(range_buf1);
		else
			f_elt++;
	}
	if (!f_is_NULL)
		return make_list_of_IRanges(range_buf2, f_levels);
	PROTECT(ans = make_CompressedIRangesList(range_buf1, OP_buf,
						 ans_breakpoints));
	UNPROTECT(2);
	return ans;
}

