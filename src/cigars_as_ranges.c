#include "cigars_as_ranges.h"

#include "IRanges_interface.h"
#include "S4Vectors_interface.h"

#include "cigar_ops_visibility.h"
#include "inspect_cigars.h"

#include <string.h>  /* for memcpy */


/* TODO: This should go in the IRanges package and be exposed
   via the IRanges C interface, like we did with
   new_list_of_IRanges_from_IntPairAEAE() */
static SEXP new_CompressedIntegerList_from_IntAEAE(const IntAEAE *aeae)
{
	size_t ans_len = IntAEAE_get_nelt(aeae);

	/* Construct 'ans_partitioning'. */
	SEXP end = PROTECT(NEW_INTEGER((R_xlen_t) ans_len));
	size_t offset = 0;
	for (size_t i = 0; i < ans_len; i++) {
		const IntAE *ae = aeae->elts[i];
		offset += IntAE_get_nelt(ae);
		INTEGER(end)[i] = offset;
	}
	SEXP ans_partitioning =
		PROTECT(new_PartitioningByEnd("PartitioningByEnd", end, NULL));

	/* Construct 'unlisted_ans'. */
	SEXP unlisted_ans = PROTECT(NEW_INTEGER((R_xlen_t) offset));
	int *unlisted_ans_p = INTEGER(unlisted_ans);
	for (size_t i = 0; i < ans_len; i++) {
		const IntAE *ae = aeae->elts[i];
		size_t ae_nelt = IntAE_get_nelt(ae);
		if (ae_nelt == 0)
			continue;
		memcpy(unlisted_ans_p, ae->elts, ae_nelt * sizeof(int));
		unlisted_ans_p += ae_nelt;
	}

	/* Construct and return 'ans'. */
	SEXP ans = PROTECT(new_CompressedList("CompressedIntegerList",
					      unlisted_ans, ans_partitioning));
	UNPROTECT(4);
	return ans;
}

static SEXP make_1col_DFrame(SEXP col, const char *colname, int nrows)
{
	/* Wrap column in named list of length 1. */
	SEXP ans_listData = PROTECT(NEW_LIST(1));
	SET_VECTOR_ELT(ans_listData, 0, col);
	SEXP ans_listData_names = PROTECT(mkString(colname));
	SET_NAMES(ans_listData, ans_listData_names);

	/* Wrap named list of length 1 in DFrame object. */
	SEXP ans_nrows = PROTECT(ScalarInteger(nrows));
	SEXP ans = PROTECT(new_DFrame("DFrame", ans_listData, ans_nrows,
				      R_NilValue));
	UNPROTECT(4);
	return ans;
}

static void set_mcol_on_IRanges(SEXP x, const char *colname, SEXP col)
{
	static SEXP elementMetadata_symbol = NULL;

	int x_len = get_IRanges_length(x);
	SEXP x_mcols = PROTECT(make_1col_DFrame(col, colname, x_len));
	if (elementMetadata_symbol == NULL)
		elementMetadata_symbol = install("elementMetadata");
	SET_SLOT(x, elementMetadata_symbol, x_mcols);
	UNPROTECT(1);
	return;
}


/****************************************************************************
 * drop_or_append_or_merge_range()
 * parse_cigar_ranges()
 */

static void drop_or_append_or_merge_range(int start, int width,
		int drop_empty_range, int merge_range, int nelt0,
		IntPairAE *range_buf,
		CharAEAE *OP_buf, char OP,
		IntAE *OPL_buf1, IntAEAE *OPL_buf2, int OPL)
{
	if (drop_empty_range && width == 0)  /* Drop. */
		return;
	int buf_nelt = IntPairAE_get_nelt(range_buf);
	if (merge_range && buf_nelt > nelt0) {
		/* The incoming range should never overlap with the previous
		   incoming range i.e. 'start' should always be > the end of
		   the previous incoming range. */
		int buf_nelt_minus_1 = buf_nelt - 1;
		int prev_end_plus_1 = range_buf->a->elts[buf_nelt_minus_1] +
				      range_buf->b->elts[buf_nelt_minus_1];
		if (start == prev_end_plus_1) {
			/* Merge. */
			range_buf->b->elts[buf_nelt_minus_1] += width;
			if (OP_buf != NULL) {
				CharAE *prev_elt =
					OP_buf->elts[buf_nelt_minus_1];
				CharAE_insert_at(prev_elt,
					CharAE_get_nelt(prev_elt), OP);
			}
			if (OPL_buf2 != NULL) {
				IntAE *prev_elt =
					OPL_buf2->elts[buf_nelt_minus_1];
				IntAE_insert_at(prev_elt,
					IntAE_get_nelt(prev_elt), OPL);
			}
			return;
		}
	}
	/* Append. */
	IntPairAE_insert_at(range_buf, buf_nelt, start, width);
	if (OP_buf != NULL) {
		CharAE *new_elt = new_CharAE(1);
		CharAE_insert_at(new_elt, 0, OP);
		CharAEAE_insert_at(OP_buf, buf_nelt, new_elt);
	}
	if (OPL_buf1 != NULL) {
		IntAE_insert_at(OPL_buf1, buf_nelt, OPL);
	}
	if (OPL_buf2 != NULL) {
		IntAE *new_elt = new_IntAE(1, 1, OPL);
		IntAEAE_insert_at(OPL_buf2, buf_nelt, new_elt);
	}
	return;
}

/* Make sure _init_ops_lkup_table() is called before parse_cigar_ranges(). */
static const char *parse_cigar_ranges(const char *cigar_string,
		int space, int lmmpos,
		int drop_empty_ranges, int reduce_ranges,
		IntPairAE *range_buf,
		CharAEAE *OP_buf, IntAE *OPL_buf1, IntAEAE *OPL_buf2)
{
	int buf_nelt0 = IntPairAE_get_nelt(range_buf);
	int cigar_offset = 0;
	int start = lmmpos;
	int n, OPL /* Operation Length */;
	char OP /* Operation */;
	while ((n = _next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return _get_cigar_parsing_error();
		int width = _op_is_visible(OP, space) ? OPL : 0;
		if (_is_in_ops(OP))
			drop_or_append_or_merge_range(start, width,
						      drop_empty_ranges,
						      reduce_ranges, buf_nelt0,
						      range_buf,
						      OP_buf, OP,
						      OPL_buf1, OPL_buf2, OPL);
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
	SEXP ans = PROTECT(new_list_of_IRanges_from_IntPairAEAE("IRanges",
								range_buf));
	//SEXP ans_names = PROTECT(duplicate(names));
	//SET_NAMES(ans, ans_names);
	SET_NAMES(ans, names);
	UNPROTECT(2);
	return ans;
}

static SEXP make_CompressedIRangesList(const IntPairAE *range_buf,
		SEXP breakpoints,
		const CharAEAE *OP_buf,
		const IntAE *OPL_buf1, const IntAEAE *OPL_buf2)
{
	SEXP unlisted_ans =
		PROTECT(new_IRanges_from_IntPairAE("IRanges", range_buf));
	if (OP_buf != NULL) {
		SEXP nms = PROTECT(new_CHARACTER_from_CharAEAE(OP_buf));
		set_IRanges_names(unlisted_ans, nms);
		UNPROTECT(1);
	}
	if (OPL_buf1 != NULL) {
		SEXP oplen = PROTECT(new_INTEGER_from_IntAE(OPL_buf1));
		set_mcol_on_IRanges(unlisted_ans, "oplen", oplen);
		UNPROTECT(1);
	} else if (OPL_buf2 != NULL) {
		SEXP oplen =
		    PROTECT(new_CompressedIntegerList_from_IntAEAE(OPL_buf2));
		set_mcol_on_IRanges(unlisted_ans, "oplen", oplen);
		UNPROTECT(1);
	}
	SEXP ans_partitioning =
		PROTECT(new_PartitioningByEnd("PartitioningByEnd",
					      breakpoints, NULL));
	SEXP ans =
		PROTECT(new_CompressedList("CompressedIRangesList",
					   unlisted_ans, ans_partitioning));
	UNPROTECT(3);
	return ans;
}

/* --- .Call ENTRY POINT ---
   Args:
     cigar:  character vector containing extended CIGAR strings.
     space:  single integer indicating one of the 8 supported spaces (defined
             at the top of the cigar_extent.c file).
     flag:   NULL or an integer vector of the same length as 'cigar'
             containing the SAM flag for each read. Serves only as a way to
             indicate whether a read is mapped or not. According to the SAM
             Spec v1.4, flag bit 0x4 is the only reliable place to tell
             whether a segment (or read) is mapped (bit is 0) or not (bit is 1).
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
SEXP C_cigars_as_ranges(SEXP cigar, SEXP space, SEXP flag, SEXP lmmpos, SEXP f,
		SEXP ops, SEXP drop_empty_ranges, SEXP reduce_ranges,
		SEXP with_ops, SEXP with_oplens)
{
	int cigar_len = LENGTH(cigar);
	const int *flag_p;
	if (flag != R_NilValue)
		flag_p = INTEGER(flag);
	_init_ops_lkup_table(ops);
	int space0 = INTEGER(space)[0];
	int pos_len = LENGTH(lmmpos);
	const int *lmmpos_p = INTEGER(lmmpos);
	int f_is_NULL = f == R_NilValue;

	IntPairAE *range_buf1;
	IntPairAEAE *range_buf2;
	SEXP ans_breakpoints;
	int *breakpoints;
	SEXP f_levels;
	const int *f_p;
	if (f_is_NULL) {
		int ans_len = cigar_len;
		/* We will typically generate at least 'cigar_len' ranges. */
		range_buf1 = new_IntPairAE(ans_len, 0);
		ans_breakpoints = PROTECT(NEW_INTEGER(ans_len));
		breakpoints = INTEGER(ans_breakpoints);
	} else {
		f_levels = GET_LEVELS(f);
		int ans_len = LENGTH(f_levels);
		range_buf2 = new_IntPairAEAE(ans_len, ans_len);
		f_p = INTEGER(f);
	}

	int drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	int reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	CharAEAE *OP_buf = NULL;
	IntAE *OPL_buf1 = NULL;
	IntAEAE *OPL_buf2 = NULL;
	if (f_is_NULL) {
		if (LOGICAL(with_ops)[0])
			OP_buf = new_CharAEAE(cigar_len, 0);
		if (LOGICAL(with_oplens)[0]) {
			if (reduce_ranges0)
				OPL_buf2 = new_IntAEAE(cigar_len, 0);
			else
				OPL_buf1 = new_IntAE(cigar_len, 0, 0);
		}
	}
	for (int i = 0; i < cigar_len; i++) {
		if (flag != R_NilValue) {
			if (*flag_p == NA_INTEGER) {
				if (f_is_NULL)
					UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (*flag_p & 0x004) {
				/* The CIGAR of an unmapped read doesn't
				   produce any range i.e. it's treated as an
				   empty CIGAR. */
				goto for_tail;
			}
		}
		SEXP cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		const char *cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		if (*lmmpos_p == NA_INTEGER || *lmmpos_p == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'lmmpos[%d]' is NA or 0", i + 1);
		}
		if (!f_is_NULL) {
			if (*f_p == NA_INTEGER)
				error("'f[%d]' is NA", i + 1);
			range_buf1 = range_buf2->elts[*f_p - 1];
		}
		const char *errmsg = parse_cigar_ranges(
					cigar_string, space0, *lmmpos_p,
					drop_empty_ranges0, reduce_ranges0,
					range_buf1, OP_buf, OPL_buf1, OPL_buf2);
		if (errmsg != NULL) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flag != R_NilValue)
			flag_p++;
		if (pos_len != 1)
			lmmpos_p++;
		if (f_is_NULL) {
			*(breakpoints++) = IntPairAE_get_nelt(range_buf1);
		} else {
			f_p++;
		}
	}
	if (!f_is_NULL)
		return make_list_of_IRanges(range_buf2, f_levels);
	SEXP ans =
		PROTECT(make_CompressedIRangesList(range_buf1, ans_breakpoints,
						   OP_buf, OPL_buf1, OPL_buf2));
	UNPROTECT(2);
	return ans;
}

