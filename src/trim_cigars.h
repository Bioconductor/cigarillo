#ifndef _TRIM_CIGARS_H_
#define _TRIM_CIGARS_H_

#include <Rdefines.h>

SEXP C_trim_cigars_along_ref(
	SEXP cigar,
	SEXP Lnpos,
	SEXP Rnpos
);

SEXP C_trim_cigars_along_query(
	SEXP cigar,
	SEXP Lnpos,
	SEXP Rnpos
);

#endif  /* _TRIM_CIGARS_H_ */

