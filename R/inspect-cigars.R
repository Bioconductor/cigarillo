### =========================================================================
### Inspect CIGAR strings
### -------------------------------------------------------------------------


validate_cigars <- function(cigars)
{
    cigars <- normarg_cigars(cigars)
    cigarillo.Call("C_validate_cigars", cigars, 0L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transform CIGARs into other useful representations
###

explode_cigar_ops <- function(cigars, ops=CIGAR_OPS)
{
    cigars <- normarg_cigars(cigars)
    ops <- normarg_ops(ops)
    cigarillo.Call("C_explode_cigar_ops", cigars, ops)
}

explode_cigar_oplens <- function(cigars, ops=CIGAR_OPS)
{
    cigars <- normarg_cigars(cigars)
    ops <- normarg_ops(ops)
    cigarillo.Call("C_explode_cigar_oplens", cigars, ops)
}

cigars_as_RleList <- function(cigars)
{
    cigar_ops <- explode_cigar_ops(cigars)
    cigar_op_lengths <- explode_cigar_oplens(cigars)
    if (length(cigars) == 0L) {
        unlisted_cigar_ops <- character(0)
        unlisted_cigar_op_lengths <- integer(0)
    } else {
        unlisted_cigar_ops <- unlist(cigar_ops, use.names=FALSE)
        unlisted_cigar_op_lengths <- unlist(cigar_op_lengths, use.names=FALSE)
    }

    ## Prepare 'ans_flesh'.
    ans_flesh <- Rle(unlisted_cigar_ops, unlisted_cigar_op_lengths)

    ## Prepare 'ans_skeleton'.
    nops_per_cigar <- elementNROWS(cigar_op_lengths)
    ans_breakpoints <- cumsum(unlisted_cigar_op_lengths)[cumsum(nops_per_cigar)]
    ans_skeleton <- PartitioningByEnd(ans_breakpoints)

    ## Relist.
    relist(ans_flesh, ans_skeleton)
}

