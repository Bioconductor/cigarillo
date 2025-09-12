#include <R_ext/Rdynload.h>

#include "inspect_cigars.h"
#include "cigar_extent.h"
#include "trim_cigars.h"
#include "cigars_as_ranges.h"
#include "position_mapping.h"
#include "map_ref_ranges_to_query.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* inspect_cigars.c */
	CALLMETHOD_DEF(C_validate_cigars, 2),
	CALLMETHOD_DEF(C_explode_cigar_ops, 2),
	CALLMETHOD_DEF(C_explode_cigar_oplens, 2),
	CALLMETHOD_DEF(C_tabulate_cigar_ops, 1),

/* cigar_extent.c */
	CALLMETHOD_DEF(C_cigar_extent, 3),

/* trim_cigars.c */
	CALLMETHOD_DEF(C_trim_cigars_along_ref, 3),
	CALLMETHOD_DEF(C_trim_cigars_along_query, 3),

/* cigars_as_ranges.c */
	CALLMETHOD_DEF(C_cigars_as_ranges, 9),

/* position_mapping.c */
	CALLMETHOD_DEF(C_query_pos_as_ref_pos, 4),
	CALLMETHOD_DEF(C_ref_pos_as_query_pos, 4),

/* map_ref_ranges_to_query.c */
	CALLMETHOD_DEF(C_map_ref_ranges_to_query, 4),

	{NULL, NULL, 0}
};

void R_init_cigarillo(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, 0);
	return;
}

