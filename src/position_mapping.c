#include "position_mapping.h"

#include "inspect_cigars.h"


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
 *   cigar:       character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns an integer vector of positions along the reference space. This
 * assumes that the query positions actually occur in the reference space,
 * outside of any deletions or insertions.
 */
SEXP C_query_pos_as_ref_pos(SEXP query_pos, SEXP cigar, SEXP lmmpos,
			    SEXP narrow_left)
{
	int npos, i;
	SEXP ref_pos;

	npos = LENGTH(query_pos);
	PROTECT(ref_pos = allocVector(INTSXP, npos));
	for (i = 0; i < npos; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigar, i));
		INTEGER(ref_pos)[i] = _to_ref(INTEGER(query_pos)[i],
					      cig_i, INTEGER(lmmpos)[i],
					      asLogical(narrow_left));
	}

	UNPROTECT(1);
	return ref_pos;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   ref_pos:     positions along the reference space
 *   cigar:       character string containing the extended CIGAR
 *   lmmpos:      1-based leftmost mapping POSition
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns an integer vector of positions along the query space. This
 * assumes that the reference positions actually occur in the query space,
 * outside of any deletions or insertions.
 */
SEXP C_ref_pos_as_query_pos(SEXP ref_pos, SEXP cigar, SEXP lmmpos,
			    SEXP narrow_left)
{
	int npos, i;
	SEXP query_pos;

	npos = LENGTH(ref_pos);
	PROTECT(query_pos = allocVector(INTSXP, npos));
	for (i = 0; i < npos; i++) {
		const char *cig_i = CHAR(STRING_ELT(cigar, i));
		INTEGER(query_pos)[i] = _to_query(INTEGER(ref_pos)[i],
						  cig_i, INTEGER(lmmpos)[i],
						  asLogical(narrow_left));
	}

	UNPROTECT(1);
	return query_pos;
}

