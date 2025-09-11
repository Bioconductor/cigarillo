### =========================================================================
### Coordinate mapping
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### query_locs_to_ref_locs()
### ref_locs_to_query_locs()
###

.normarg_query_locs <- function(query_locs)
{
    if (!is.numeric(query_locs))
        stop(wmsg("'query_locs' must be a vector of integers"))
    if (!is.integer(query_locs))
        query_locs <- as.integer(query_locs)
    query_locs
}

.normarg_ref_locs <- function(ref_locs)
{
    if (!is.numeric(ref_locs))
        stop(wmsg("'ref_locs' must be a vector of integers"))
    if (!is.integer(ref_locs))
        ref_locs <- as.integer(ref_locs)
    ref_locs
}

### Returns an integer vector parallel to 'query_locs'.
query_locs_to_ref_locs <- function(query_locs, cigars, lmmpos, narrow.left)
{
    query_locs <- .normarg_query_locs(query_locs)
    cigars <- normarg_cigars(cigars)
    if (length(query_locs) != length(cigars))
        stop(wmsg("'query_locs' and 'cigars' must have the same length"))
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    if (!isTRUEorFALSE(narrow.left))
        stop(wmsg("'narrow.left' must be TRUE or FALSE"))
    cigarillo.Call("C_query_locs_to_ref_locs",
                   query_locs, cigars, lmmpos, narrow.left)
}

### Returns an integer vector parallel to 'ref_locs'.
ref_locs_to_query_locs <- function(ref_locs, cigars, lmmpos, narrow.left)
{
    ref_locs <- .normarg_ref_locs(ref_locs)
    cigars <- normarg_cigars(cigars)
    if (length(ref_locs) != length(cigars))
        stop(wmsg("'ref_locs' and 'cigars' must have the same length"))
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    if (!isTRUEorFALSE(narrow.left))
        stop(wmsg("'narrow.left' must be TRUE or FALSE"))
    cigarillo.Call("C_ref_locs_to_query_locs",
                   ref_locs, cigars, lmmpos, narrow.left)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### queryLoc2refLoc()
### queryLocs2refLocs()
###
### What were those supposed to do again? Aren't they redundant with
### query_locs_to_ref_locs()?
###

### Not exported.
queryLoc2refLoc <- function(qloc, cigars, lmmpos=1L)
{
    stop(wmsg("not implemented yet"))
}

### Not exported.
queryLocs2refLocs <- function(qlocs, cigars, lmmpos=1L, flags=NULL)
{
    stop(wmsg("not implemented yet"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### map_query_locs_to_ref_locs()
### map_ref_locs_to_query_locs()
###
### Both functions return a list of length 4:
###   1. start of local query position
###   2. end of local query position
###   3. index of 'start' used in match ('from_hits')
###   4. index of 'lmmpos' used in match ('to_hits')
### All list elements are integer vectors. This assumes that the reference
### positions actually occur in the read alignment region, outside of
### any deletions or insertions.
###

map_query_locs_to_ref_locs <- function(start, end, cigars, lmmpos)
{
    cigars <- normarg_cigars(cigars)
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    cigarillo.Call("C_map_query_locs_to_ref_locs", start, end, cigars, lmmpos)
}

map_ref_locs_to_query_locs <- function(start, end, cigars, lmmpos)
{
    cigars <- normarg_cigars(cigars)
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    cigarillo.Call("C_map_ref_locs_to_query_locs", start, end, cigars, lmmpos)
}

