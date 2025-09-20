#ifndef _CIGARS_AS_RANGES_H_
#define _CIGARS_AS_RANGES_H_

#include <Rdefines.h>

SEXP C_cigars_as_ranges(
	SEXP cigars,
	SEXP space,
	SEXP flags,
	SEXP lmmpos,
	SEXP f,
	SEXP ops,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges,
	SEXP with_ops,
	SEXP with_oplens
);

#endif  /* _CIGARS_AS_RANGES_H_ */

