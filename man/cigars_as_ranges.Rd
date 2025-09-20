\name{cigars_as_ranges}

\alias{cigars_as_ranges}

\alias{cigars_as_ranges_along_ref}
\alias{cigars_as_ranges_along_query}
\alias{cigars_as_ranges_along_pwa}

\title{Turn CIGAR strings into ranges of positions}

\description{
  Turn CIGAR strings into ranges of positions relative to
  the "query space", "reference space", or "pairwise alignment space".
}

\usage{
cigars_as_ranges_along_ref(cigars,
          N.regions.removed=FALSE,
          flags=NULL, lmmpos=1L, f=NULL,
          ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
          with.ops=FALSE, with.oplens=FALSE)

cigars_as_ranges_along_query(cigars,
          before.hard.clipping=FALSE, after.soft.clipping=FALSE,
          flags=NULL,
          ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
          with.ops=FALSE, with.oplens=FALSE)

cigars_as_ranges_along_pwa(cigars,
          N.regions.removed=FALSE, dense=FALSE,
          flags=NULL,
          ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
          with.ops=FALSE, with.oplens=FALSE)
}

\arguments{
  \item{cigars}{
    A character vector (or factor) containing CIGAR strings.
  }
  \item{N.regions.removed}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigars_as_ranges_along_ref} reports ranges
    with respect to the "reference space" from which the N regions have
    been removed, and \code{cigars_as_ranges_along_pwa} reports them with
    respect to the "pairwise alignment space" from which the N regions
    have been removed.
  }
  \item{flags}{
    \code{NULL} or an integer vector parallel to \code{cigars} that
    contains the SAM/BAM flags corresponding to each CIGAR string.

    According to the SAM Spec v1.4, flag bit 0x4 is the only reliable place
    to tell whether a segment (or read) is mapped (bit is 0) or not (bit
    is 1). If the \code{flags} argument is supplied, then
    \code{cigars_as_ranges_along_ref}, \code{cigars_as_ranges_along_query},
    and \code{cigars_as_ranges_along_pwa} don't produce any range
    for unmapped reads i.e. they treat them as if their CIGAR was empty
    (independently of what their CIGAR is).
  }
  \item{lmmpos}{
    An integer vector containing the 1-based leftmost mapping POSition
    of each alignment with respect to the "reference space". These are
    the 1-based leftmost positions/coordinates of each (eventually clipped)
    query sequence with respect to the subject.

    \code{lmmpos} must be a single integer, or an integer vector parallel
    to \code{cigars}.
  }
  \item{f}{
    \code{NULL} or a factor parallel to \code{cigars}.

    If \code{NULL} (the default), then the ranges are grouped by alignment
    i.e. the returned \link[IRanges]{IRangesList} object has 1 list element
    per element in \code{cigars}. Otherwise they are grouped by factor level
    i.e. the returned \link[IRanges]{IRangesList} object has 1 list element
    per level in \code{f} and is named with those levels.

    For example, if \code{f} is a factor containing the chromosome for each
    read, then the returned \link[IRanges]{IRangesList} object will have
    1 list element per chromosome and each list element will contain all
    the ranges on that chromosome.
  }
  \item{ops}{
    Character vector where the elements are single letters representing
    valid CIGAR operations. Must be a subset of \code{\link{CIGAR_OPS}}.
    See \code{?\link{CIGAR_OPS}} for more information.

    Only the operations listed in \code{ops} will be turned into ranges.
  }
  \item{drop.empty.ranges}{
    \code{TRUE} or \code{FALSE}.

    Should empty ranges be dropped?
  }
  \item{reduce.ranges}{
    \code{TRUE} or \code{FALSE}.

    Should adjacent ranges coming from the same cigar be merged or not?
    Using \code{TRUE} can significantly reduce the size of the returned
    object.
  }
  \item{with.ops}{
    \code{TRUE} or \code{FALSE}.

    Should the returned ranges be named/labeled with their corresponding
    CIGAR operation? Only supported when \code{f} is \code{NULL}.
  }
  \item{with.oplens}{
    \code{TRUE} or \code{FALSE}.

    If \code{with.oplens} is \code{TRUE}, then the returned
    \link[IRanges]{IRangesList} object will carry the lengths of the
    CIGAR operations in an inner metadata column named \code{oplen}.
    Only supported when \code{f} is \code{NULL}.
  }
  \item{before.hard.clipping}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigars_as_ranges_along_query} reports
    ranges with respect to the "query space" to which the H regions
    have been added.
    Note that \code{before.hard.clipping} and \code{after.soft.clipping}
    cannot both be \code{TRUE}.
  }
  \item{after.soft.clipping}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigars_as_ranges_along_query} reports
    ranges with respect to the "query space" from which the S regions
    have been removed.
    Note that \code{before.hard.clipping} and \code{after.soft.clipping}
    cannot both be \code{TRUE}.
  }
  \item{dense}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigars_as_ranges_along_pwa} reports
    ranges with respect to the "pairwise alignment space" from which
    the I, D, and N regions have been removed.
    Note that \code{N.regions.removed} and \code{dense} cannot both
    be \code{TRUE}.
  }
}

\value{
  An \link[IRanges]{IRangesList} object (more precisely a
  \link[IRanges]{CompressedIRangesList} object) with one list
  element per element in \code{cigars}.

  However, if \code{f} is a factor, then the returned
  \link[IRanges]{IRangesList} object returned by
  \code{cigars_as_ranges_along_ref()} is a
  \link[IRanges]{SimpleIRangesList} object (instead of
  \link[IRanges]{CompressedIRangesList}). In that case it has one
  list element per level in \code{f}, and is named with those levels.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{cigar_ops_visibility}} for an introduction to CIGAR
          operations and their visibility in various "projection spaces".

    \item \link{inspect_cigars} for basic utility functions to parse
          and inspect CIGAR strings.

    \item \code{\link{tabulate_cigar_ops}} to count the occurences of CIGAR
           operations in a vector of CIGAR strings.

    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \code{\link{trim_cigars_along_ref}} and
          \code{\link{trim_cigars_along_query}} to trim CIGAR strings
          along the "reference space" and "query space", respectively.

    \item \code{\link{project_positions}} to project positions from query
          to reference space and vice versa.

    \item \code{\link{project_sequences}} to project sequences from one
          space to the other.

    \item The \link[IRanges]{IRanges} and \link[IRanges]{IRangesList}
          classes in the \pkg{IRanges} package.
  }
}

\examples{
cigar1 <- "3H15M55N4M2I6M2D5M6S"
my_cigars <- c("40M2I9M", cigar1, "2S10M2000N15M", "3H33M5H")

## ---------------------------------------------------------------------
## Turn CIGAR strings into ranges along the "reference space"
## ---------------------------------------------------------------------

cigars_as_ranges_along_ref(cigar1, with.ops=TRUE, with.oplens=TRUE)[[1]]

cigars_as_ranges_along_ref(cigar1, reduce.ranges=TRUE,
                           with.ops=TRUE, with.oplens=TRUE)[[1]]

ops <- setdiff(CIGAR_OPS, "N")

cigars_as_ranges_along_ref(cigar1, ops=ops,
                           with.ops=TRUE, with.oplens=TRUE)[[1]]

cigars_as_ranges_along_ref(cigar1, ops=ops, reduce.ranges=TRUE,
                           with.ops=TRUE, with.oplens=TRUE)[[1]]

ops <- setdiff(CIGAR_OPS, c("D", "N"))

cigars_as_ranges_along_ref(cigar1, ops=ops,
                           with.ops=TRUE, with.oplens=TRUE)[[1]]

lmmpos <- c(1, 1001, 1,  351)

cigars_as_ranges_along_ref(my_cigars, lmmpos=lmmpos,
                           with.ops=TRUE, with.oplens=TRUE)

cigars_as_ranges_along_ref(my_cigars, lmmpos=lmmpos,
                           ops=setdiff(CIGAR_OPS, "N"),
                           reduce.ranges=TRUE)

cigars_as_ranges_along_ref(my_cigars, lmmpos=lmmpos,
                           ops=setdiff(CIGAR_OPS, c("D", "N")),
                           reduce.ranges=TRUE)

seqnames <- factor(c("chr6", "chr6", "chr2", "chr6"),
                   levels=c("chr2", "chr6"))
ops <- c("M", "=", "X", "I", "D")
cigars_as_ranges_along_ref(my_cigars, lmmpos=lmmpos, f=seqnames, ops=ops)

## ---------------------------------------------------------------------
## Turn CIGAR strings into ranges along the "query space"
## ---------------------------------------------------------------------

cigars_as_ranges_along_query(my_cigars, with.ops=TRUE, with.oplens=TRUE)

## ---------------------------------------------------------------------
## Turn CIGAR strings into ranges along the "pairwise alignment space"
## ---------------------------------------------------------------------

cigars_as_ranges_along_pwa(my_cigars, with.ops=TRUE, with.oplens=TRUE)
cigars_as_ranges_along_pwa(my_cigars, dense=TRUE,
                           with.ops=TRUE, with.oplens=TRUE)
}

\keyword{manip}
