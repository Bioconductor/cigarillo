### =========================================================================
### Tabulate CIGAR operations
### -------------------------------------------------------------------------


tabulate_cigar_ops <- function(cigars, oplens.as.weights=FALSE)
{
    cigars <- normarg_cigars(cigars)
    if (!isTRUEorFALSE(oplens.as.weights))
        stop(wmsg("'oplens.as.weights' must be TRUE or FALSE"))
    ans <- cigarillo.Call("C_tabulate_cigar_ops", cigars, oplens.as.weights)
    stopifnot(identical(CIGAR_OPS, colnames(ans)))  # sanity check
    ans
}

