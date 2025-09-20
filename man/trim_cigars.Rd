\name{trim_cigars}

\alias{trim_cigars}

\alias{trim_cigars_along_ref}
\alias{trim_cigars_along_query}
\alias{narrow_cigars_along_ref}
\alias{narrow_cigars_along_query}

\title{Trim CIGAR strings along reference or query space}

\description{
  The CIGAR string associated with a pairwise alignment describes
  the alignment in its entirety in the sense that it covers all the
  positions in the alignment. However, there might be situations
  where one is only interested in a particular portion of the alignment,
  that is, in the portion of the alignment that is left after trimming
  it by a given number of positions on its left and/or right ends.
  Furthermore, one might want to know the effect of this trimming
  on the original CIGAR string.

  The \pkg{cigarillo} package provides two core functions,
  \code{trim_cigars_along_ref} and \code{trim_cigars_along_query},
  to compute the CIGAR string that describes a "trimmed alignment".
  Both take:
  \itemize{
    \item the original CIGAR string i.e. the CIGAR string that
          describes the alignment before trimming
    \item the numbers of left/right positions to trim
  }
  Both functions return the "trimmed CIGAR string", that is, the
  CIGAR string that describes the "trimmed alignment".

  The only difference between the two function is how the numbers of
  left and right positions to trim are counted: with respect to the
  "reference space" for \code{trim_cigars_along_ref}, and with respect
  to the "query space" for \code{trim_cigars_along_query}.

  Both functions are vectorized.
}

\usage{
trim_cigars_along_ref(cigars, Lnpos=0L, Rnpos=0L)
trim_cigars_along_query(cigars, Lnpos=0L, Rnpos=0L)

## Wrappers to the above that do the same thing but via
## the "narrow()" interface:
narrow_cigars_along_ref(cigars, start=NA, end=NA, width=NA)
narrow_cigars_along_query(cigars, start=NA, end=NA, width=NA)
}

\arguments{
  \item{cigars}{
    A character vector (or factor) containing CIGAR strings.
  }
  \item{Lnpos,Rnpos}{
    The numbers of left/right positions to trim.

    Each of \code{Lnpos} and \code{Rnpos} must be a non-negative integer,
    or a vector of non-negative integers of the same length as \code{cigars}.

    Note that the numbers of left and right positions to trim are counted
    with respect to the "reference space" for \code{trim_cigars_along_ref},
    and with respect to the "query space" for \code{trim_cigars_along_query}.
  }
  \item{start,end,width}{
    Vectors of integers. \code{NA}s and negative values are allowed
    and "solved" similarly to what \code{IRanges::narrow()} does.
    See \code{?IRanges::\link[IRanges]{narrow}} in the \pkg{IRanges}
    package for more information.
  }
}

\value{
  A character vector of the same length as \code{cigars} that contains
  the "trimmed CIGAR strings".

  In addition the vector has an "rshift" attribute which is an integer
  vector of the same length as \code{cigars}. It contains the values that
  would need to be added to the POS field of a SAM/BAM file as a
  consequence of this trimming.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{cigar_ops_visibility}} for an introduction to CIGAR
          operations and their visibility in various "projection spaces".

    \item \code{\link{tabulate_cigar_ops}} to count the occurences of CIGAR
           operations in a vector of CIGAR strings.

    \item \link{inspect_cigars} for basic utility functions to parse
          and inspect CIGAR strings.

    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \link{cigars_as_ranges} to turn CIGAR strings into ranges
          of positions.

    \item \code{\link{project_positions}} to project positions from query
          to reference space and vice versa.

    \item \code{\link{project_sequences}} to project sequences from one
          space to the other.
  }
}

\examples{
cigar1 <- "3H15M55N4M2I6M2D5M6S"

## trim_cigars_along_ref():
trim_cigars_along_ref(cigar1)  # only drops the soft/hard clipping
trim_cigars_along_ref(cigar1, Lnpos=9)
trim_cigars_along_ref(cigar1, Lnpos=14)
trim_cigars_along_ref(cigar1, Lnpos=14, Rnpos=16)
trim_cigars_along_ref(cigar1, Lnpos=15)
#trim_cigars_along_ref(cigar1, Lnpos=15, Rnpos=17)  # error! (empty cigar)
trim_cigars_along_ref(cigar1, Lnpos=70)
trim_cigars_along_ref(cigar1, Lnpos=71)
trim_cigars_along_ref(cigar1, Lnpos=74)

## trim_cigars_along_query():
trim_cigars_along_query(cigar1, Lnpos=3, Rnpos=2)
trim_cigars_along_query(cigar1, Lnpos=9)
trim_cigars_along_query(cigar1, Lnpos=18)
trim_cigars_along_query(cigar1, Lnpos=23)

## Using the "narrow()" interface:

stopifnot(
  ## narrow_cigars_along_ref() vs trim_cigars_along_ref():
  identical(narrow_cigars_along_ref(cigar1, start=10),
            trim_cigars_along_ref(cigar1, Lnpos=9)),
  identical(narrow_cigars_along_ref(cigar1, start=15),
            trim_cigars_along_ref(cigar1, Lnpos=14)),
  identical(narrow_cigars_along_ref(cigar1, start=15, width=57),
            trim_cigars_along_ref(cigar1, Lnpos=14, Rnpos=16)),
  identical(narrow_cigars_along_ref(cigar1, start=16),
            trim_cigars_along_ref(cigar1, Lnpos=15)),
  identical(narrow_cigars_along_ref(cigar1, start=71),
            trim_cigars_along_ref(cigar1, Lnpos=70)),
  identical(narrow_cigars_along_ref(cigar1, start=72),
            trim_cigars_along_ref(cigar1, Lnpos=71)),
  identical(narrow_cigars_along_ref(cigar1, start=75),
            trim_cigars_along_ref(cigar1, Lnpos=74)),

  ## narrow_cigars_along_query() vs trim_cigars_along_query():
  identical(narrow_cigars_along_query(cigar1, start=4, end=-3),
            trim_cigars_along_query(cigar1, Lnpos=3, Rnpos=2)),
  identical(narrow_cigars_along_query(cigar1, start=10),
            trim_cigars_along_query(cigar1, Lnpos=9)),
  identical(narrow_cigars_along_query(cigar1, start=19),
            trim_cigars_along_query(cigar1, Lnpos=18)),
  identical(narrow_cigars_along_query(cigar1, start=24),
            trim_cigars_along_query(cigar1, Lnpos=23))
)
}

\keyword{manip}
