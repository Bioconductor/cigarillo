### =========================================================================
### map_ref_ranges_to_query()
### -------------------------------------------------------------------------
###
### Highly specialized utility functions whose main purpose is to
### support GenomicAlignments:::.mapFromAlignments().
### Too specialized to be of general interest!
###


.normarg_start <- function(start, what="start")
{
    if (!is.numeric(start))
        stop(wmsg("'", what, "' must be an integer vector"))
    if (!is.integer(start))
        start <- as.integer(start)
    if (anyNA(start))
        stop(wmsg("'", what, "' cannot contain NAs"))
    start
}

### start, end:    two parallel integer vectors describing ranges along the
###                reference space (input ranges);
### cigar, lmmpos: two parallel vectors (one character, one integer).
###
### Finds the hits between the input ranges and the vector of
### cigar/lmmpos pairs. An input range is considered to have a hit with
### a cigar/lmmpos pair if it's a valid range within the extent of the
### corresponding alignment along the reference space.
###
### Returns the hits in a 4-column data.frame with 1 hit per row.
### Note that the rows are sorted by "from_hit" first then by "to_hit".
map_ref_ranges_to_query <- function(start, end, cigars, lmmpos)
{
    start <- .normarg_start(start)
    end <- .normarg_start(end, what="end")
    if (length(start) != length(end))
        stop(wmsg("'start' and 'end' must have the same length"))
    cigars <- normarg_cigars(cigars)
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    C_ans <- cigarillo.Call("C_map_ref_ranges_to_query",
                            start, end, cigars, lmmpos)
    structure(C_ans, names=c("start", "end", "from_hit", "to_hit"),
              class="data.frame", row.names=seq_along(C_ans[[1L]]))
}

### A reimplementation of map_ref_ranges_to_query() that is based on
### findOverlaps(). Hundreds times faster than map_ref_ranges_to_query()
### for medium size input (i.e. when nb of input ranges x nb of cigars is
### between 1e6 and 250e6). Thousands to hundreds of thousands times faster
### or more for big inputs (i.e. when nb of input ranges x nb of cigars
### is > 500e6).
###
### Note that rows in the returned data.frame are only guaranteed to be
### sorted by "from_hit". Use 'strictly.sort.hits=TRUE' to have them
### sorted by "from_hit" first then by "to_hit". Note that in this case,
### the returned data.frame is identical to the data.frame returned by
### map_ref_ranges_to_query().
fast_map_ref_ranges_to_query <- function(start, end, cigars, lmmpos,
                                         strictly.sort.hits=FALSE)
{
    input_ranges <- IRanges(start, end)
    cigar_ranges <- IRanges(lmmpos, width=cigar_extent_along_ref(cigars))
    if (!isTRUEorFALSE(strictly.sort.hits))
        stop(wmsg("'strictly.sort.hits' must be TRUE or FALSE"))
    hits <- findOverlaps(input_ranges, cigar_ranges, type="within")
    ## findOverlaps() always returns a SortedByQueryHits object so the hits
    ## are already guaranteed to be sorted by query hit. However,
    ## if 'strictly.sort.hits' is TRUE, we want them sorted by query hit
    ## first then by subject hit. Quick and easy way to achieve this is
    ## to transpose the Hits object twice. This is very fast!
    if (strictly.sort.hits)
        hits <- t(t(hits))
    ans_from_hit <- queryHits(hits)
    ans_to_hit <- subjectHits(hits)
    ref_ranges <- input_ranges[ans_from_hit]
    cig <- cigars[ans_to_hit]
    pos <- lmmpos[ans_to_hit]
    ans_start <- ref_pos_as_query_pos(start(ref_ranges), cig, pos, FALSE)
    ans_end <- ref_pos_as_query_pos(end(ref_ranges), cig, pos, TRUE)
    data.frame(
        start=ans_start,
        end=ans_end,
        from_hit=ans_from_hit,
        to_hit=ans_to_hit)
}

