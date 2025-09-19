### =========================================================================
### Turn CIGAR strings into ranges of positions
### -------------------------------------------------------------------------


cigars_as_ranges <-
    function(cigars, space,
             flags=NULL, lmmpos=1L, f=NULL,
             ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
             with.ops=FALSE, with.oplens=FALSE)
{
    cigars <- normarg_cigars(cigars)
    if (!isSingleNumber(space))
        stop("'space' must be a single integer")
    if (!is.integer(space))
        space <- as.integer(space)
    flags <- normarg_flags(flags, cigars)
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    if (!is.null(f)) {
        if (!is.factor(f))
            stop("'f' must be NULL or a factor")
        if (length(f) != length(cigars))
            stop("'f' must have the same length as 'cigars'")
    }
    ops <- normarg_ops(ops)
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(reduce.ranges))
        stop("'reduce.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.ops))
        stop("'with.ops' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.oplens))
        stop("'with.oplens' must be TRUE or FALSE")
    cigarillo.Call("C_cigars_as_ranges",
                   cigars, space, flags, lmmpos, f,
                   ops, drop.empty.ranges, reduce.ranges,
                   with.ops, with.oplens)
}

cigars_as_ranges_along_ref <-
    function(cigars, N.regions.removed=FALSE,
             flags=NULL, lmmpos=1L, f=NULL,
             ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
             with.ops=FALSE, with.oplens=FALSE)
{
    space <- select_reference_space(N.regions.removed)
    C_ans <- cigars_as_ranges(cigars, space, flags, lmmpos, f,
                              ops, drop.empty.ranges, reduce.ranges,
                              with.ops, with.oplens)
    if (is.null(f))
        return(C_ans)
    compress <- length(C_ans) >= 200L
    IRangesList(C_ans, compress=compress)
}

cigars_as_ranges_along_query <-
    function(cigars, before.hard.clipping=FALSE, after.soft.clipping=FALSE,
             flags=NULL,
             ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
             with.ops=FALSE, with.oplens=FALSE)
{
    space <- select_query_space(before.hard.clipping, after.soft.clipping)
    cigars_as_ranges(cigars, space, flags, 1L, NULL,
                     ops, drop.empty.ranges, reduce.ranges,
                     with.ops, with.oplens)
}

cigars_as_ranges_along_pwa <-
    function(cigars, N.regions.removed=FALSE, dense=FALSE,
             flags=NULL,
             ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
             with.ops=FALSE, with.oplens=FALSE)
{
    space <- select_pairwise_space(N.regions.removed, dense)
    cigars_as_ranges(cigars, space, flags, 1L, NULL,
                     ops, drop.empty.ranges, reduce.ranges,
                     with.ops, with.oplens)
}

