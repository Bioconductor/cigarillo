\name{cigar_extent}

\alias{cigar_extent}

\alias{cigar_extent_along_ref}
\alias{cigar_extent_along_query}
\alias{cigar_extent_along_pwa}

\title{
  Calculate the number of positions spanned by a CIGAR
}

\description{
  The \emph{extent} (or length) of a pairwise alignment is the number of
  positions that it spans. Note that positions can be counted with respect
  to the "reference space", "query space", or "pairwise alignment space".
  This means that the \emph{extent} of a pairwise alignment depends on the
  space that we're using to count positions.

  The \emph{extent} of a CIGAR string is simply the \emph{extent} of the
  alignment that it describes.

  The \pkg{cigarillo} package provides three functions to calculate the
  \emph{extent} of a CIGAR string:
  \itemize{
    \item \code{cigar_extent_along_ref} calculates the extent along
          the "reference space".
    \item \code{cigar_extent_along_query} calculates the extent along
          the "query space".
    \item \code{cigar_extent_along_pwa} calculates the extent along
          the "pairwise alignment space".
  }

  The three functions are vectorized.
}

\usage{
cigar_extent_along_ref(cigars, flags=NULL,
        N.regions.removed=FALSE)

cigar_extent_along_query(cigars, flags=NULL,
        before.hard.clipping=FALSE, after.soft.clipping=FALSE)

cigar_extent_along_pwa(cigars, flags=NULL,
        N.regions.removed=FALSE, dense=FALSE)
}

\arguments{
  \item{cigars}{
    A character vector or factor containing CIGAR strings.
  }
  \item{flags}{
    \code{NULL} or an integer vector containing the SAM flag for
    each read.

    According to the SAM Spec v1.4, flag bit 0x4 is the only reliable place
    to tell whether a segment (or read) is mapped (bit is 0) or not (bit
    is 1). If the \code{flags} argument is supplied, then
    \code{cigar_extent_along_ref}, \code{cigar_extent_along_query}, and
    \code{cigar_extent_along_pwa} return \code{NA}s for unmapped reads.
  }
  \item{N.regions.removed}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigar_extent_along_ref} reports the CIGAR
    extents with respect to the "reference space" from which the N regions
    have been removed, and \code{cigar_extent_along_pwa} reports them with
    respect to the "pairwise alignment space" from which the N regions
    have been removed.
  }
  \item{before.hard.clipping}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigar_extent_along_query} reports the CIGAR
    extents with respect to the "query space" to which the H regions
    have been added.
    Note that \code{before.hard.clipping} and \code{after.soft.clipping}
    cannot both be \code{TRUE}.
  }
  \item{after.soft.clipping}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigar_extent_along_query} reports the CIGAR
    extents with respect to the "query space" from which the S regions
    have been removed.
    Note that \code{before.hard.clipping} and \code{after.soft.clipping}
    cannot both be \code{TRUE}.
  }
  \item{dense}{
    \code{TRUE} or \code{FALSE}.

    If \code{TRUE}, then \code{cigar_extent_along_pwa} reports the CIGAR
    extents with respect to the "pairwise alignment space" from which
    the I, D, and N regions have been removed.
    Note that \code{N.regions.removed} and \code{dense} cannot both
    be \code{TRUE}.
  }
}

\value{
  For \code{cigar_extent_along_ref} and \code{cigar_extent_along_pwa}:
  An integer vector of the same length as \code{cigars} where each
  element is the extent of the alignment with respect to the reference
  and pairwise space, respectively.
  More precisely, for \code{cigar_extent_along_ref}, the returned
  extents are the lengths of the alignments on the reference,
  N gaps included (except if \code{N.regions.removed} is \code{TRUE}).
  NAs or \code{"*"} in \code{cigars} will produce NAs in the returned vector.

  For \code{cigar_extent_along_query}: An integer vector of the same
  length as \code{cigars} where each element is the length of the
  corresponding query sequence as inferred from the CIGAR string.
  Note that, by default (i.e. if \code{before.hard.clipping} and
  \code{after.soft.clipping} are \code{FALSE}), this is the length
  of the query sequence stored in the SAM/BAM file.
  If \code{before.hard.clipping} or \code{after.soft.clipping}
  is \code{TRUE}, the returned extents are the lengths of the query
  sequences before hard clipping or after soft clipping.
  NAs or \code{"*"} in \code{cigars} will produce NAs in the returned vector.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{inspect_cigars} for basic utility functions to parse
          and inspect CIGAR strings.

    \item \code{\link{trim_cigars_along_ref}} and
          \code{\link{trim_cigars_along_query}} to trim CIGAR strings
          along the "reference space" and "query space", respectively.

    \item \link{cigars_as_ranges} for turning CIGAR strings into ranges
          of positions.
  }
}

\examples{
my_cigars <- c("40M2I9M", "3H15M55N4M2I6M2D5M6S",
               "2S10M2000N15M", "3H33M5H")

## Extents along the "reference space":
cigar_extent_along_ref(my_cigars)

## Extents along the "query space":
cigar_extent_along_query(my_cigars)
cigar_extent_along_query(my_cigars, before.hard.clipping=TRUE)

## Extents along the "pairwise alignment space":
cigar_extent_along_pwa(my_cigars)
cigar_extent_along_pwa(my_cigars, dense=TRUE)
}

\keyword{manip}
