#ifndef _COORDINATE_MAPPING_H_
#define _COORDINATE_MAPPING_H_

#include <Rdefines.h>

SEXP C_query_locs_to_ref_locs(
	SEXP query_locs,
	SEXP cigar,
	SEXP lmmpos,
	SEXP narrow_left
);

SEXP C_ref_locs_to_query_locs(
	SEXP ref_locs,
	SEXP cigar,
	SEXP lmmpos,
	SEXP narrow_left
);

SEXP C_map_query_locs_to_ref_locs(
	SEXP start,
	SEXP end,
	SEXP cigar,
	SEXP lmmpos
);

SEXP C_map_ref_locs_to_query_locs(
	SEXP start,
	SEXP end,
	SEXP cigar,
	SEXP lmmpos
);

#endif  /* _COORDINATE_MAPPING_H_ */

