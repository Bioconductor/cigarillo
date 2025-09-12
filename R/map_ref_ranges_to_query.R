### =========================================================================
### map_ref_ranges_to_query()
### -------------------------------------------------------------------------
###
### This is a highly specialized utility whose only purpose is to
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
    structure(C_ans, names=c("start", "end", "range_idx", "cigar_idx"),
              class="data.frame", row.names=seq_along(C_ans[[1L]]))
}

