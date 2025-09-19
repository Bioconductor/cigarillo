### =========================================================================
### Project positions from query to reference space, and vice versa
### -------------------------------------------------------------------------
###


### Does the same job as normarg_lmmpos() really. Unify?
.normarg_pos <- function(pos, cigars, what)
{
    if (!is.numeric(pos))
        stop(wmsg("'", what, "' must be an integer vector"))
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (anyNA(pos))
        stop(wmsg("'", what, "' cannot contain NAs"))
    if (length(pos) != length(cigars))
        stop(wmsg("'", what, "' and 'cigars' must have the same length"))
    pos
}

### Returns an integer vector parallel to 'query_pos'.
query_pos_as_ref_pos <- function(query_pos, cigars, lmmpos, narrow.left)
{
    cigars <- normarg_cigars(cigars)
    query_pos <- .normarg_pos(query_pos, cigars, "query_pos")
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    if (!isTRUEorFALSE(narrow.left))
        stop(wmsg("'narrow.left' must be TRUE or FALSE"))
    cigarillo.Call("C_query_pos_as_ref_pos",
                   query_pos, cigars, lmmpos, narrow.left)
}

### Returns an integer vector parallel to 'ref_pos'.
ref_pos_as_query_pos <- function(ref_pos, cigars, lmmpos, narrow.left)
{
    cigars <- normarg_cigars(cigars)
    ref_pos <- .normarg_pos(ref_pos, cigars, "ref_pos")
    lmmpos <- normarg_lmmpos(lmmpos, cigars)
    if (!isTRUEorFALSE(narrow.left))
        stop(wmsg("'narrow.left' must be TRUE or FALSE"))
    cigarillo.Call("C_ref_pos_as_query_pos",
                   ref_pos, cigars, lmmpos, narrow.left)
}

