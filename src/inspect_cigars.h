#ifndef _INSPECT_CIGARS_H_
#define _INSPECT_CIGARS_H_

#include <Rdefines.h>

const char *_get_cigar_parsing_error();

int _next_cigar_OP(
	const char *cigar_string,
	int offset,
	char *OP,
	int *OPL
);

int _prev_cigar_OP(
	const char *cigar_string,
	int offset,
	char *OP,
	int *OPL
);

void _init_ops_lkup_table(SEXP ops);

int _is_in_ops(char OP);

SEXP C_validate_cigars(
	SEXP cigars,
	SEXP ans_type
);

SEXP C_explode_cigar_ops(
	SEXP cigars,
	SEXP ops
);

SEXP C_explode_cigar_oplens(
	SEXP cigars,
	SEXP ops
);

#endif  /* _INSPECT_CIGARS_H_ */

