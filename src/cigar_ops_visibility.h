#ifndef _CIGAR_OPS_VISIBILITY_H_
#define _CIGAR_OPS_VISIBILITY_H_

#include <Rdefines.h>

int _op_is_visible(
	char OP,
	int space
);

SEXP C_cigar_ops_visibility(SEXP ops);

#endif  /* _CIGAR_OPS_VISIBILITY_H_ */

