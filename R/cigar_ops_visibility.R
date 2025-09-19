### =========================================================================
### Visibility of CIGAR operations
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cigar_ops_visibility()
###

### The 8 "projection spaces" below are also defined at the top of the
### src/cigar_extent.c file.
PROJECTION_SPACES <- c(
    "reference",
    "reference-N-regions-removed",
    "query",
    "query-before-hard-clipping",
    "query-after-soft-clipping",
    "pairwise",
    "pairwise-N-regions-removed",
    "pairwise-dense"
)

cigar_ops_visibility <- function(ops=CIGAR_OPS)
{
    ops <- normarg_ops(ops)
    ans <- cigarillo.Call("C_cigar_ops_visibility", ops)
    dimnames(ans) <- list(PROJECTION_SPACES, ops)
    ans
}

