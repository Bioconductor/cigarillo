### =========================================================================
### Some low-level internal utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


normarg_cigars <- function(cigars)
{
    if (is.factor(cigars))
        cigars <- as.character(cigars)
    if (!is.character(cigars))
        stop(wmsg("'cigars' must be a character vector or factor"))
    cigars
}

### lmmpos: 1-based leftmost mapping POSition
### See https://samtools.github.io/hts-specs/SAMv1.pdf
normarg_lmmpos <- function(lmmpos, cigars)
{
    if (!is.numeric(lmmpos))
        stop(wmsg("'lmmpos' must be a vector of integers"))
    if (!is.integer(lmmpos))
        lmmpos <- as.integer(lmmpos)
    if (length(lmmpos) != 1L && length(lmmpos) != length(cigars))
        stop(wmsg("'lmmpos' must have length 1 or ",
                  "the same length as 'cigars'"))
    lmmpos
}

cigarillo.Call <- function(.NAME, ...) .Call2(.NAME, ..., PACKAGE="cigarillo")

