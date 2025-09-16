#include "map_ref_ranges_to_query.h"

#include "S4Vectors_interface.h"

#include "project_positions.h"


/*
 * Code in this file originally written by Valerie Obenchain in Dec 2014
 * for the GenomicAlignments package.
 * Code copied from GenomicAlignments to cigarillo on Sep 12, 2025.
 *
 * Note that the current implementation uses nested for loops to find all
 * the hits between the input ranges and cigar/lmmpos pairs, which is very
 * inefficient.
 * See fast_map_ref_ranges_to_query() in R/map_ref_ranges_to_query.R for a
 * more efficient approach based on findOverlaps().
 */


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   start, end:    two parallel integer vectors describing ranges along
 *                  the reference space (input ranges);
 *   cigar, lmmpos: two parallel vectors (one character, one integer).
 * Returns a list of length four that describes the hits between the input
 * ranges and the cigar/lmmpos pairs. All list elements are parallel integer
 * vectors of length N, where N is the number of hits.
 * The four list elements are:
 *   1. start of reference range relative to query space
 *   2. end of reference range relative to query space
 *   3. index of input range involved in hit
 *   4. index of cigar/lmmpos pair involved in hit
 * Note that an input range is considered to have a hit with a cigar/lmmpos
 * pair if it's a valid range within the extent of the corresponding
 * alignment along the reference space.
 */
SEXP C_map_ref_ranges_to_query(SEXP start, SEXP end, SEXP cigar, SEXP lmmpos)
{
	SEXP ans, ans_start, ans_end, ans_qhits, ans_shits;
	IntAE *sbuf, *ebuf, *qhbuf, *shbuf;
	int i, j, s, e;
	int nranges = LENGTH(start);
	int ncigar = LENGTH(cigar);

	sbuf = new_IntAE(0, 0, 0);
	ebuf = new_IntAE(0, 0, 0);
	qhbuf = new_IntAE(0, 0, 0);
	shbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < nranges; i++) {
		for (j = 0; j < ncigar; j++) {
			const char *cig_j = CHAR(STRING_ELT(cigar, j));
			int pos_j = INTEGER(lmmpos)[j];
			s = _to_query(INTEGER(start)[i], cig_j, pos_j, FALSE);
			if (s == NA_INTEGER)
				continue;
			e = _to_query(INTEGER(end)[i], cig_j, pos_j, TRUE);
			if (e == NA_INTEGER)
				continue;
			IntAE_insert_at(sbuf, IntAE_get_nelt(sbuf), s);
			IntAE_insert_at(ebuf, IntAE_get_nelt(ebuf), e);
			IntAE_insert_at(qhbuf, IntAE_get_nelt(qhbuf), i + 1);
			IntAE_insert_at(shbuf, IntAE_get_nelt(shbuf), j + 1);
		}
	}

	PROTECT(ans = NEW_LIST(4));
	PROTECT(ans_start = new_INTEGER_from_IntAE(sbuf));
	PROTECT(ans_end = new_INTEGER_from_IntAE(ebuf));
	PROTECT(ans_qhits = new_INTEGER_from_IntAE(qhbuf));
	PROTECT(ans_shits = new_INTEGER_from_IntAE(shbuf));
	SET_VECTOR_ELT(ans, 0, ans_start);
	SET_VECTOR_ELT(ans, 1, ans_end);
	SET_VECTOR_ELT(ans, 2, ans_qhits);
	SET_VECTOR_ELT(ans, 3, ans_shits);
	UNPROTECT(5);
	return ans;
}

