#ifndef _TABULATE_CIGAR_OPS_H_
#define _TABULATE_CIGAR_OPS_H_

#include <Rdefines.h>

SEXP C_tabulate_cigar_ops(
	SEXP cigars,
	SEXP oplens_as_weights
);

#endif  /* _TABULATE_CIGAR_OPS_H_ */

