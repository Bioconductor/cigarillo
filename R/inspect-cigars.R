### =========================================================================
### Inspect CIGAR strings
### -------------------------------------------------------------------------


### See p. 4 of the SAM Spec v1.4 at http://samtools.sourceforge.net/ for the
### list of CIGAR operations and their meanings.
CIGAR_OPS <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")

normarg_ops <- function(ops)
{
    if (is.null(ops))
        return(ops)
    if (!is.character(ops))
        stop("'ops' must be a character vector")
    if (any(is.na(ops)))
        stop("'ops' cannot contain NAs")
    if (length(ops) == 1L) {
        ops <- strsplit(ops, NULL, fixed=TRUE)[[1L]]
    } else if (any(nchar(ops) != 1L)) {
        stop("when 'length(ops) != 1', all its elements ",
             "must be single letters")
    }
    if (anyDuplicated(ops))
        stop("'ops' cannot contain duplicated letters")
    if (!all(ops %in% CIGAR_OPS))
        stop("'ops' contains invalid CIGAR operations")
    ops
}

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tabulate_cigar_ops()
###

tabulate_cigar_ops <- function(cigars)
{
    cigars <- normarg_cigars(cigars)
    ans <- cigarillo.Call("C_tabulate_cigar_ops", cigars)
    stopifnot(identical(CIGAR_OPS, colnames(ans)))  # sanity check
    ans
}

