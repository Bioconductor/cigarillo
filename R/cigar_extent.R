### =========================================================================
### Extent of a CIGAR = nb of positions it spans
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers used by the cigar_extent_along_<space>() and
### cigars_as_ranges_along_<space>() functions
###

select_reference_space <- function(N.regions.removed)
{
    if (!isTRUEorFALSE(N.regions.removed))
        stop(wmsg("'N.regions.removed' must be TRUE or FALSE"))
    if (N.regions.removed) {
        space <- 2L  # REFERENCE_N_REGIONS_REMOVED
    } else {
        space <- 1L  # REFERENCE
    }
    space
}

select_query_space <- function(before.hard.clipping, after.soft.clipping)
{
    if (!isTRUEorFALSE(before.hard.clipping))
        stop(wmsg("'before.hard.clipping' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(after.soft.clipping))
        stop(wmsg("'after.soft.clipping' must be TRUE or FALSE"))
    if (before.hard.clipping) {
        if (after.soft.clipping)
            stop(wmsg("'before.hard.clipping' and 'after.soft.clipping' ",
                      "cannot both be TRUE"))
        space <- 4L  # QUERY_BEFORE_HARD_CLIPPING
    } else if (after.soft.clipping) {
        space <- 5L  # QUERY_AFTER_SOFT_CLIPPING
    } else {
        space <- 3L  # QUERY
    }
    space
}

select_pairwise_space <- function(N.regions.removed, dense)
{
    if (!isTRUEorFALSE(N.regions.removed))
        stop(wmsg("'N.regions.removed' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(dense))
        stop(wmsg("'dense' must be TRUE or FALSE"))
    if (N.regions.removed) {
        if (dense)
           stop(wmsg("'N.regions.removed' and 'dense' ",
                     "cannot both be TRUE"))
        space <- 7L  # PAIRWISE_N_REGIONS_REMOVED
    } else if (dense) {
        space <- 8L  # PAIRWISE_DENSE
    } else {
        space <- 6L  # PAIRWISE
    }
    space
}

normarg_flags <- function(flags, cigars)
{
    if (!is.null(flags)) {
        if (!is.numeric(flags))
            stop(wmsg("'flags' must be NULL or a vector of integers"))
        if (!is.integer(flags))
            flags <- as.integer(flags)
        if (length(cigars) != length(flags))
            stop(wmsg("'cigars' and 'flags' must have the same length"))
    }
    flags
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The cigar_extent_along_<space>() functions
###

.cigar_extent <- function(cigars, space, flags)
{
    cigars <- normarg_cigars(cigars)
    flags <- normarg_flags(flags, cigars)
    if (!isSingleNumber(space))
        stop(wmsg("'space' must be a single integer"))
    if (!is.integer(space))
        space <- as.integer(space)
    cigarillo.Call("C_cigar_extent", cigars, space, flags)
}

cigar_extent_along_ref <- function(cigars,
                                   N.regions.removed=FALSE,
                                   flags=NULL)
{
    space <- select_reference_space(N.regions.removed)
    .cigar_extent(cigars, space, flags)
}

cigar_extent_along_query <- function(cigars,
                                     before.hard.clipping=FALSE,
                                     after.soft.clipping=FALSE,
                                     flags=NULL)
{
    space <- select_query_space(before.hard.clipping, after.soft.clipping)
    .cigar_extent(cigars, space, flags)
}

cigar_extent_along_pwa <- function(cigars,
                                   N.regions.removed=FALSE, dense=FALSE,
                                   flags=NULL)
{
    space <- select_pairwise_space(N.regions.removed, dense)
    .cigar_extent(cigars, space, flags)
}

