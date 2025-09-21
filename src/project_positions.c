#include "project_positions.h"

#include "explode_cigars.h"


/*
 * Code in this file originally written by Michael Lawrence in 2012 for
 * the GenomicRanges package.
 * Code moved from GenomicRanges to GenomicAlignments on Dec 6, 2013.
 * Code copied from GenomicAlignments to cigarillo on Sep 12, 2025.
 */


/* Turns single position along query space ('query_pos') into position
   along reference space.
   If 'query_pos' cannot be mapped NA is returned. */
int _to_ref(int query_pos, const char *cig, int lmmpos, Rboolean narrow_left)
{
  int ref_pos = query_pos + lmmpos - 1;
  int n, offset = 0, OPL, query_consumed = 0;
  char OP;

  while (query_consumed < query_pos &&
         (n = _next_cigar_OP(cig, offset, &OP, &OPL)))
  {
    switch (OP) {
      /* Alignment match (can be a sequence match or mismatch) */
      case 'M': case '=': case 'X':
          query_consumed += OPL;
          break;
      /* Insertion to the reference */
      case 'I': {
        int width_from_insertion_start = query_pos - query_consumed;
        Rboolean query_pos_past_insertion = width_from_insertion_start > OPL;
        if (query_pos_past_insertion) {
          ref_pos -= OPL;
        } else {
          ref_pos -= width_from_insertion_start;
          if (!narrow_left) {
            ref_pos += 1;
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
        ref_pos += OPL;
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
    ref_pos = NA_INTEGER;

  return ref_pos;
}

/* Turns single position along reference space ('ref_pos') into position
   along query space.
   If 'ref_pos' cannot be mapped NA is returned. */
int _to_query(int ref_pos, const char *cig, int lmmpos, Rboolean narrow_left)
{
  int query_pos = ref_pos - lmmpos + 1;
  int n, offset = 0, OPL, query_consumed = 0;
  char OP;

  while (query_consumed < query_pos &&
         (n = _next_cigar_OP(cig, offset, &OP, &OPL)))
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
      query_pos += OPL;
      query_consumed += OPL;
      break;
    /* Deletion from the reference */
    case 'D':
    /* Skipped region from reference; narrow to query */
    case 'N':
      {
        Rboolean query_pos_past_gap = query_pos - query_consumed > OPL;
        if (query_pos_past_gap) {
          query_pos -= OPL;
        } else {
          if (narrow_left) {
            query_pos = query_consumed;
          } else {
            query_pos = query_consumed + 1;
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

  if (query_pos <= 0 || n == 0)
    query_pos = NA_INTEGER;

  return query_pos;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   query_pos:   positions along the query space
 *   cigars:      character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns an integer vector of positions along the reference space. This
 * assumes that the query positions actually occur in the reference space,
 * outside of any deletions or insertions.
 */
SEXP C_query_pos_as_ref_pos(SEXP query_pos, SEXP cigars, SEXP lmmpos,
			    SEXP narrow_left)
{
	int npos = LENGTH(query_pos);
	SEXP ref_pos = PROTECT(allocVector(INTSXP, npos));
	int lmmpos_len = LENGTH(lmmpos);
	const int *lmmpos_p = INTEGER(lmmpos);
	for (int i = 0; i < npos; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigars, i));
		INTEGER(ref_pos)[i] = _to_ref(INTEGER(query_pos)[i],
					      cig_i, *lmmpos_p,
					      asLogical(narrow_left));
		if (lmmpos_len != 1)
			lmmpos_p++;
	}

	UNPROTECT(1);
	return ref_pos;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   ref_pos:     positions along the reference space
 *   cigars:      character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns an integer vector of positions along the query space. This
 * assumes that the reference positions actually occur in the query space,
 * outside of any deletions or insertions.
 */
SEXP C_ref_pos_as_query_pos(SEXP ref_pos, SEXP cigars, SEXP lmmpos,
			    SEXP narrow_left)
{
	int npos = LENGTH(ref_pos);
	SEXP query_pos = PROTECT(allocVector(INTSXP, npos));
	int lmmpos_len = LENGTH(lmmpos);
	const int *lmmpos_p = INTEGER(lmmpos);
	for (int i = 0; i < npos; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigars, i));
		INTEGER(query_pos)[i] = _to_query(INTEGER(ref_pos)[i],
						  cig_i, *lmmpos_p,
						  asLogical(narrow_left));
		if (lmmpos_len != 1)
			lmmpos_p++;
	}

	UNPROTECT(1);
	return query_pos;
}

