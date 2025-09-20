### =========================================================================
### Project sequences from one space to the other
### -------------------------------------------------------------------------
###
### This is a complete rewrite of GenomicAlignments::sequenceLayer().
### Better, simpler implementation, that is also about twice faster!
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_indel_ops()
###

.get_indel_ops <- function(from, to)
{
    stopifnot(isSingleInteger(from), isSingleInteger(to))
    all_indel_ops <- c("I", "D", "N", "S", "H")
    vis_mat <- cigar_ops_visibility(all_indel_ops)
    all_indel_ops[vis_mat[from, ] != vis_mat[to, ]]
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

### Returns an XStringSetList derivative parallel to IRangesList
### object 'indel_at'.
.make_replacement_value <- function(indel_at, class,
                                    I.letter, D.letter, N.letter,
                                    S.letter, H.letter)
{
    stopifnot(is(indel_at, "CompressedIRangesList"), isSingleString(class))
    unlisted_indel_at <- unlist(indel_at, use.names=FALSE)
    unlisted_indel_at_names <- names(unlisted_indel_at)
    stopifnot(!is.null(unlisted_indel_at_names))

    unlisted_ans <- rep.int(as("", class), length(unlisted_indel_at))
    names(unlisted_ans) <- unlisted_indel_at_names

    inject_idx <- which(width(unlisted_indel_at) == 0L)
    if (length(inject_idx) != 0L) {
        oplens <- mcols(unlisted_indel_at)[inject_idx, "oplen"]
        ops <- unlisted_indel_at_names[inject_idx]
        filler_seqs <- .make_filler_sequences(oplens, ops,
                                              seqtype(unlisted_ans),
                                              I.letter, D.letter, N.letter,
                                              S.letter, H.letter)
        stopifnot(identical(names(filler_seqs), ops))
        unlisted_ans[inject_idx] <- filler_seqs
    }
    relist(unlisted_ans, indel_at)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### project_sequences()
###

project_sequences <- function(x, cigars, from="query", to="reference",
                              I.letter="-", D.letter="-", N.letter=".",
                              S.letter="+", H.letter="+")
{
    if (!is(x, "XStringSet"))
        stop(wmsg("'x' must be an XStringSet object"))
    cigars <- normarg_cigars(cigars)
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
        stop(wmsg("'H.letter' must be the same as 'S.letter' ",
                  "when 'from' is not \"query\" and 'to' ",
                  "is \"query-before-hard-clipping\""))

    from <- match(from, PROJECTION_SPACES)
    to <- match(to, PROJECTION_SPACES)
    indel_ops <- .get_indel_ops(from, to)
    indel_at <- cigars_as_ranges(cigars, from, ops=indel_ops,
                                 with.ops=TRUE, with.oplens=TRUE)
    value <- .make_replacement_value(indel_at, class(x),
                                     I.letter, D.letter, N.letter,
                                     S.letter, H.letter)
    replaceAt(x, indel_at, value=value)
}

