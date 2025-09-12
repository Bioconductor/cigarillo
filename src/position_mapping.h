#ifndef _POSITION_MAPPING_H_
#define _POSITION_MAPPING_H_

#include <Rdefines.h>

int _to_ref(
	int query_pos,
	const char *cig,
	int lmmpos,
	Rboolean narrow_left
);

int _to_query(
	int ref_pos,
	const char *cig,
	int lmmpos,
	Rboolean narrow_left
);

SEXP C_query_pos_as_ref_pos(
	SEXP query_pos,
	SEXP cigar,
	SEXP lmmpos,
	SEXP narrow_left
);

SEXP C_ref_pos_as_query_pos(
	SEXP ref_pos,
	SEXP cigar,
	SEXP lmmpos,
	SEXP narrow_left
);

#endif  /* _POSITION_MAPPING_H_ */

