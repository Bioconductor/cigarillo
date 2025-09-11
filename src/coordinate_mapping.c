#include "coordinate_mapping.h"

#include "S4Vectors_interface.h"

#include "inspect_cigars.h"


/* Returns integer position of 'query_loc' mapped to genome-based space.
   If 'query_loc' cannot be mapped NA is returned. */
static int to_ref(int query_loc, const char *cig0, int lmmpos,
		  Rboolean narrow_left)
{
  int ref_loc = query_loc + lmmpos - 1;
  int n, offset = 0, OPL, query_consumed = 0;
  char OP;

  while (query_consumed < query_loc &&
         (n = _next_cigar_OP(cig0, offset, &OP, &OPL)))
  {
    switch (OP) {
      /* Alignment match (can be a sequence match or mismatch) */
      case 'M': case '=': case 'X':
          query_consumed += OPL;
          break;
      /* Insertion to the reference */
      case 'I': {
        int width_from_insertion_start = query_loc - query_consumed;
        Rboolean query_loc_past_insertion = width_from_insertion_start > OPL;
        if (query_loc_past_insertion) {
          ref_loc -= OPL;
        } else {
          ref_loc -= width_from_insertion_start;
          if (!narrow_left) {
            ref_loc += 1;
          }
        }
        query_consumed += OPL;
        break;
      }
      /* Soft clip on the read */
      case 'S':
        query_consumed += OPL;
        break;
      /* Deletion from the reference */
      case 'D':
      case 'N': /* Skipped region from reference; narrow to query */
        ref_loc += OPL;
        break;
      /* Hard clip on the read */
      case 'H':
        break;
      /* Silent deletion from the padded reference */
      case 'P':
        break;
      default:
        break;
    }
    offset += n;
  }

  if (n == 0)
    ref_loc = NA_INTEGER;

  return ref_loc;
}

/* Returns integer position of 'ref_loc' mapped to local space.
   If 'ref_loc' cannot be mapped NA is returned. */
static int to_query(int ref_loc, const char *cig0, int lmmpos,
		    Rboolean narrow_left)
{
  int query_loc = ref_loc - lmmpos + 1;
  int n, offset = 0, OPL, query_consumed = 0;
  char OP;

  while (query_consumed < query_loc &&
         (n = _next_cigar_OP(cig0, offset, &OP, &OPL)))
  {
    switch (OP) {
    /* Alignment match (can be a sequence match or mismatch) */
    case 'M': case '=': case 'X':
      query_consumed += OPL;
      break;
    /* Insertion to the reference */
    case 'I':
    /* Soft clip on the read */
    case 'S':
      query_loc += OPL;
      query_consumed += OPL;
      break;
    /* Deletion from the reference */
    case 'D':
    /* Skipped region from reference; narrow to query */
    case 'N':
      {
        Rboolean query_loc_past_gap = query_loc - query_consumed > OPL;
        if (query_loc_past_gap) {
          query_loc -= OPL;
        } else {
          if (narrow_left) {
            query_loc = query_consumed;
          } else {
            query_loc = query_consumed + 1;
          }
        }
      }
      break;
    /* Hard clip on the read */
    case 'H':
      break;
    /* Silent deletion from the padded reference */
    case 'P':
      break;
    default:
      break;
    }
    offset += n;
  }

  if (query_loc <= 0 || n == 0)
    query_loc = NA_INTEGER;

  return query_loc;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   query_locs:  local positions in the read that we will map
 *   cigar:       character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns an integer vector of local query positions. This assumes
 * that the reference positions actually occur in the read alignment region,
 * outside of any deletions or insertions.
 */
SEXP C_query_locs_to_ref_locs(SEXP query_locs, SEXP cigar, SEXP lmmpos,
			      SEXP narrow_left)
{
	int nlocs, i;
	SEXP ref_locs;

	nlocs = LENGTH(query_locs);
	PROTECT(ref_locs = allocVector(INTSXP, nlocs));
	for (i = 0; i < nlocs; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigar, i));
		INTEGER(ref_locs)[i] = to_ref(INTEGER(query_locs)[i],
					      cig_i, INTEGER(lmmpos)[i],
					      asLogical(narrow_left));
	}

	UNPROTECT(1);
	return ref_locs;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   ref_locs:    global positions in the reference to map
 *   cigar:       character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 *
 * Returns an integer vector of local query positions. This assumes the
 * reference positions actually occur in the read alignment region,
 * outside of any deletions or insertions.
 */
SEXP C_ref_locs_to_query_locs(SEXP ref_locs, SEXP cigar, SEXP lmmpos,
			      SEXP narrow_left)
{
	int nlocs, i;
	SEXP query_locs;

	nlocs = LENGTH(ref_locs);
	PROTECT(query_locs = allocVector(INTSXP, nlocs));
	for (i = 0; i < nlocs; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigar, i));
		INTEGER(query_locs)[i] = to_query(INTEGER(ref_locs)[i],
						  cig_i, INTEGER(lmmpos)[i],
						  asLogical(narrow_left));
	}

	UNPROTECT(1);
	return query_locs;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Returns a list of length 4:
 *      - start of local query position
 *      - end of local query position
 *      - index of 'start' used in match ('from_hits')
 *      - index of 'lmmpos' used in match ('to_hits')
 * All list elements are integer vectors. This assumes that the reference
 * positions actually occur in the read alignment region, outside of
 * any deletions or insertions.
 */
SEXP C_map_query_locs_to_ref_locs(SEXP start, SEXP end, SEXP cigar, SEXP lmmpos)
{
	SEXP ans, ans_start, ans_end, ans_qhits, ans_shits;
	IntAE *sbuf, *ebuf, *qhbuf, *shbuf;
	int i, j, s, e, nlocs, ncigar;

	nlocs = LENGTH(start);
	ncigar = LENGTH(cigar);
	sbuf = new_IntAE(0, 0, 0);
	ebuf = new_IntAE(0, 0, 0);
	qhbuf = new_IntAE(0, 0, 0);
	shbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < nlocs; i++) {
		for (j = 0; j < ncigar; j++) {
			const char *cig_j = CHAR(STRING_ELT(cigar, j));
			int pos_j = INTEGER(lmmpos)[j];
			s = to_ref(INTEGER(start)[i], cig_j, pos_j, FALSE);
			if (s == NA_INTEGER)
				break;
			e = to_ref(INTEGER(end)[i], cig_j, pos_j, TRUE);
			if (e == NA_INTEGER)
				break;
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


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Returns a list of length 4:
 *   - start of local query position
 *   - end of local query position
 *   - index of 'start' used in match ('from_hits')
 *   - index of 'lmmpos' used in match ('to_hits')
 * All list elements are integer vectors. This assumes that the reference
 * positions actually occur in the read alignment region, outside of
 * any deletions or insertions.
 */
SEXP C_map_ref_locs_to_query_locs(SEXP start, SEXP end, SEXP cigar, SEXP lmmpos)
{
	SEXP ans, ans_start, ans_end, ans_qhits, ans_shits;
	IntAE *sbuf, *ebuf, *qhbuf, *shbuf;
	int i, j, s, e;
	int nlocs = LENGTH(start);
	int ncigar = LENGTH(cigar);

	sbuf = new_IntAE(0, 0, 0);
	ebuf = new_IntAE(0, 0, 0);
	qhbuf = new_IntAE(0, 0, 0);
	shbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < nlocs; i++) {
		for (j = 0; j < ncigar; j++) {
			const char *cig_j = CHAR(STRING_ELT(cigar, j));
			int pos_j = INTEGER(lmmpos)[j];
			s = to_query(INTEGER(start)[i], cig_j, pos_j, FALSE);
			if (s == NA_INTEGER)
				continue;
			e = to_query(INTEGER(end)[i], cig_j, pos_j, TRUE);
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

