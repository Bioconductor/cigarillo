### =========================================================================
### Trim CIGAR strings
### -------------------------------------------------------------------------


.normarg_npos <- function(npos, cigars, what="Lnpos")
{
    if (!is.numeric(npos))
        stop(wmsg("'", what, "' must be an integer vector"))
    if (!is.integer(npos))
        npos <- as.integer(npos)
    if (length(npos) != length(cigars)) {
        if (length(npos) != 1L)
            stop(wmsg("'", what, "' must have length 1 or ",
                      "the same length as 'cigars'"))
        npos <- rep.int(npos, length(cigars))
    }
    if (anyNA(npos))
        stop(wmsg("'", what, "' cannot contain NAs"))
    if (suppressWarnings(min(npos)) < 0L)
        stop(wmsg("'", what, "' cannot contain negative values"))
    npos
}

trim_cigars_along_ref <- function(cigars, Lnpos=0L, Rnpos=0L)
{
    cigars <- normarg_cigars(cigars)
    Lnpos <- .normarg_npos(Lnpos, cigars, what="Lnpos")
    Rnpos <- .normarg_npos(Rnpos, cigars, what="Rnpos")
    C_ans <- cigarillo.Call("C_trim_cigars_along_ref", cigars, Lnpos, Rnpos)
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}

trim_cigars_along_query <- function(cigars, Lnpos=0L, Rnpos=0L)
{
    cigars <- normarg_cigars(cigars)
    Lnpos <- .normarg_npos(Lnpos, cigars, what="Lnpos")
    Rnpos <- .normarg_npos(Rnpos, cigars, what="Rnpos")
    C_ans <- cigarillo.Call("C_trim_cigars_along_query", cigars, Lnpos, Rnpos)
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}


### Wrappers that do the same thing but via the "narrow()" interface.

.compute_LRnpos <- function(extents, start=NA, end=NA, width=NA)
{
    extents_as_ranges <- IRanges(1L, extents)
    threeranges <- threebands(extents_as_ranges,
                              start=start, end=end, width=width)
    list(Lnpos=width(threeranges$left), Rnpos=width(threeranges$right))
}

narrow_cigars_along_ref <- function(cigars, start=NA, end=NA, width=NA)
{
    extents <- cigar_extent_along_ref(cigars)
    LRnpos <- .compute_LRnpos(extents, start=start, end=end, width=width)
    trim_cigars_along_ref(cigars, LRnpos[[1L]], LRnpos[[2L]])
}

narrow_cigars_along_query <- function(cigars, start=NA, end=NA, width=NA)
{
    extents <- cigar_extent_along_query(cigars)
    LRnpos <- .compute_LRnpos(extents, start=start, end=end, width=width)
    trim_cigars_along_query(cigars, LRnpos[[1L]], LRnpos[[2L]])
}

