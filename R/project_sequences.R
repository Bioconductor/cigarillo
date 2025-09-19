### =========================================================================
### Project sequences from one space to the other
### -------------------------------------------------------------------------
###
### This is a complete rewrite of GenomicAlignments::sequenceLayer().
### Better, simpler implementation, that is also about twice faster!
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 8 .ranges_to_replace_from_<space>() helpers
###

.ranges_to_replace_from_reference <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference-N-regions-removed" = "N",
        "query"                       = c("I", "D", "N", "S"),
        "query-before-hard-clipping"  = c("I", "D", "N", "S", "H"),
        "query-after-soft-clipping"   = c("I", "D", "N"),
        "pairwise"                    = "I",
        "pairwise-N-regions-removed"  = c("I", "N"),
        "pairwise-dense"              = c("D", "N"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_ref(cigars, ops=ops,
                               with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_reference_N_regions_removed <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = "N",
        "query"                       = c("I", "D", "S"),
        "query-before-hard-clipping"  = c("I", "D", "S", "H"),
        "query-after-soft-clipping"   = c("I", "D"),
        "pairwise"                    = c("I", "N"),
        "pairwise-N-regions-removed"  = "I",
        "pairwise-dense"              = "D",
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_ref(cigars, N.regions.removed=TRUE, ops=ops,
                               with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_query <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = c("I", "D", "N", "S"),
        "reference-N-regions-removed" = c("I", "D", "S"),
        "query-before-hard-clipping"  = "H",
        "query-after-soft-clipping"   = "S",
        "pairwise"                    = c("D", "N", "S"),
        "pairwise-N-regions-removed"  = c("D", "S"),
        "pairwise-dense"              = c("I", "S"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_query(cigars, ops=ops,
                                 with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_query_before_hard_clipping <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = c("I", "D", "N", "S", "H"),
        "reference-N-regions-removed" = c("I", "D", "S", "H"),
        "query"                       = "H",
        "query-after-soft-clipping"   = c("S", "H"),
        "pairwise"                    = c("D", "N", "S", "H"),
        "pairwise-N-regions-removed"  = c("D", "S", "H"),
        "pairwise-dense"              = c("I", "S", "H"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_query(cigars, before.hard.clipping=TRUE, ops=ops,
                                 with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_query_after_soft_clipping <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = c("I", "D", "N"),
        "reference-N-regions-removed" = c("I", "D"),
        "query"                       = "S",
        "query-before-hard-clipping"  = c("S", "H"),
        "pairwise"                    = c("D", "N"),
        "pairwise-N-regions-removed"  = "D",
        "pairwise-dense"              = "I",
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_query(cigars, after.soft.clipping=TRUE, ops=ops,
                                 with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_pairwise <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = "I",
        "reference-N-regions-removed" = c("I", "N"),
        "query"                       = c("D", "N", "S"),
        "query-before-hard-clipping"  = c("D", "N", "S", "H"),
        "query-after-soft-clipping"   = c("D", "N"),
        "pairwise-N-regions-removed"  = "N",
        "pairwise-dense"              = c("I", "D", "N"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_pwa(cigars, ops=ops,
                               with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_pairwise_N_regions_removed <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = c("I", "N"),
        "reference-N-regions-removed" = "I",
        "query"                       = c("D", "S"),
        "query-before-hard-clipping"  = c("D", "S", "H"),
        "query-after-soft-clipping"   = "D",
        "pairwise"                    = "N",
        "pairwise-dense"              = c("I", "D"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_pwa(cigars, N.regions.removed=TRUE, ops=ops,
                               with.ops=TRUE, with.oplens=TRUE)
}

.ranges_to_replace_from_pairwise_dense <- function(cigars, to)
{
    stopifnot(isSingleString(to))
    ops <- switch(to,
        "reference"                   = c("D", "N"),
        "reference-N-regions-removed" = "D",
        "query"                       = c("I", "S"),
        "query-before-hard-clipping"  = c("I", "S", "H"),
        "query-after-soft-clipping"   = "I",
        "pairwise"                    = c("I", "D", "N"),
        "pairwise-N-regions-removed"  = c("I", "D"),
        stop(wmsg(to, ": invalid space")))
    cigars_as_ranges_along_pwa(cigars, dense=TRUE, ops=ops,
                               with.ops=TRUE, with.oplens=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_replacement_value()
###

### Returns a named XStringSet derivative of seqtype 'seqtype' that is
### parallel to 'oplens' and 'ops'. The names on it is the 'ops' vector.
.make_filler_sequences <- function(oplens, ops, seqtype,
                                   I.letter, D.letter, N.letter,
                                   S.letter, H.letter)
{
    stopifnot(is.integer(oplens), is.character(ops),
              length(oplens) == length(ops),
              isSingleString(seqtype))
    max_oplen_per_op <- max(splitAsList(oplens, ops))
    unique_ops <- names(max_oplen_per_op)
    biggest_fillers <- lapply(setNames(seq_along(max_oplen_per_op), unique_ops),
            function(i) {
                op <- unique_ops[[i]]
                letter <- switch(op, I=I.letter,
                                     D=D.letter,
                                     N=N.letter,
                                     S=S.letter,
                                     H=H.letter,
                                     stop(wmsg(op, ": invalid CIGAR op")))
                rep.int(letter, max_oplen_per_op[[i]])
            })
    biggest_fillers <- as(biggest_fillers, paste0(seqtype, "StringSet"))
    narrow(biggest_fillers[ops], 1L, oplens)
}

### Returns an XStringSetList derivative parallel to IRangesList object 'at'.
.make_replacement_value <- function(at, class, I.letter, D.letter, N.letter,
                                               S.letter, H.letter)
{
    stopifnot(is(at, "CompressedIRangesList"), isSingleString(class))
    unlisted_at <- unlist(at, use.names=FALSE)
    unlisted_at_names <- names(unlisted_at)
    stopifnot(!is.null(unlisted_at_names))

    unlisted_ans <- rep.int(as("", class), length(unlisted_at))
    names(unlisted_ans) <- unlisted_at_names

    inject_idx <- which(width(unlisted_at) == 0L)
    if (length(inject_idx) != 0L) {
        oplens <- mcols(unlisted_at)[inject_idx, "oplen"]
        ops <- unlisted_at_names[inject_idx]
        filler_seqs <- .make_filler_sequences(oplens, ops,
                                              seqtype(unlisted_ans),
                                              I.letter, D.letter, N.letter,
                                              S.letter, H.letter)
        stopifnot(identical(names(filler_seqs), ops))
        unlisted_ans[inject_idx] <- filler_seqs
    }
    relist(unlisted_ans, at)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### project_sequences()
###

project_sequences <- function(x, cigars, from="query", to="reference",
                              I.letter="-", D.letter="-", N.letter=".",
                              S.letter="+", H.letter="+")
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    from <- match.arg(from, PROJECTION_SPACES)
    to <- match.arg(to, PROJECTION_SPACES)
    I.letter <- Biostrings:::.normarg_padding.letter(I.letter, seqtype(x))
    D.letter <- Biostrings:::.normarg_padding.letter(D.letter, seqtype(x))
    N.letter <- Biostrings:::.normarg_padding.letter(N.letter, seqtype(x))
    S.letter <- Biostrings:::.normarg_padding.letter(S.letter, seqtype(x))
    H.letter <- Biostrings:::.normarg_padding.letter(H.letter, seqtype(x))
    if (from == to)
        return(x)

    ## Right now, the way 'S.letter' and 'H.letter' are injected in 'x' when
    ## 'to' is "query-before-hard-clipping" can result in padding in the
    ## wrong order (i.e. padding with 'H.letter' followed by padding with
    ## 'S.letter') so we temporarily work around this by enforcing 'S.letter'
    ## and 'H.letter' to be the same.
    if (from != "query" && to == "query-before-hard-clipping" &&
        as.character(S.letter) != as.character(H.letter))
        stop("'H.letter' must be the same as 'S.letter' ",
             "when 'from' is not \"query\" and 'to' ",
             "is \"query-before-hard-clipping\"")

    FUN <- paste0(".ranges_to_replace_from_", chartr("-", "_", from))
    at <- do.call(FUN, list(cigars=cigars, to=to))
    value <- .make_replacement_value(at, class(x),
                                     I.letter, D.letter, N.letter,
                                     S.letter, H.letter)
    replaceAt(x, at, value=value)
}

