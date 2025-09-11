#ifndef _CIGAR_EXTENT_H_
#define _CIGAR_EXTENT_H_

#include <Rdefines.h>

int _is_visible_in_space(
	char OP,
	int space
);

SEXP C_cigar_extent(
	SEXP cigar,
	SEXP flag,
	SEXP space
);

#endif  /* _CIGAR_EXTENT_ */

