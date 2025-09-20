\name{map_ref_ranges_to_query}

\alias{map_ref_ranges_to_query}
\alias{fast_map_ref_ranges_to_query}

\title{Map ranges relative to reference space to query space}

\description{
  Highly specialized utility functions whose main purpose is to
  support the \code{\link[GenomicAlignments]{mapToAlignments}}
  methods defined in the \pkg{GenomicAlignments} package.
  Only of interest to the authors/maintainers of these methods,
  and not really meant to be used by the end user.
}

\usage{
map_ref_ranges_to_query(start, end, cigars, lmmpos)

fast_map_ref_ranges_to_query(start, end, cigars, lmmpos,
                             strictly.sort.hits=FALSE)
}

\arguments{
  \item{start, end}{
    Two parallel integer vectors containing the starts/ends of the
    ranges to map to the "query space". Note that the positions in the
    two vectors are expected to be relative to the "reference space".
  }
  \item{cigars}{
    A character vector (or factor) containing CIGAR strings.
  }
  \item{lmmpos}{
    An integer vector parallel to \code{cigars}. For each CIGAR string
    in \code{cigars}, \code{lmmpos} must contain the 1-based leftmost
    mapping POSition of the alignment described by the CIGAR string.
    Note that these positions must be relative to the "reference space".
  }
  \item{strictly.sort.hits}{
    Whether the rows in the data.frame returned by
    \code{fast_map_ref_ranges_to_query()} should be sorted by
    \code{from_hit} first then by \code{to_hit} instead of by
    \code{from_hit} only.
    Note that when \code{strictly.sort.hits} is set to \code{TRUE},
    \code{fast_map_ref_ranges_to_query()} is guaranted to return
    the exact same data.frame as \code{map_ref_ranges_to_query()}.
  }
}

\details{
  \code{map_ref_ranges_to_query()} uses a naive and inefficient approach
  to find hits between the input ranges and the ranges implicitly defined
  by the (\code{cigars[j]}, \code{lmmpos[j]}) pairs.

  \code{fast_map_ref_ranges_to_query()} is just a reimplementation
  of \code{map_ref_ranges_to_query()} that is based on
  \code{\link[IRanges]{findOverlaps}()}. It's hundreds times faster
  than \code{map_ref_ranges_to_query()} for medium size input (i.e. when
  nb of input ranges x nb of cigars is between 1e6 and 250e6), and thousands
  to hundreds of thousands times faster or more for big inputs (i.e. when
  nb of input ranges x nb of cigars is > 500e6).
}

\value{
  A 4-column data.frame with 1 hit per row. The columns are:
  \itemize{
    \item \code{start}, \code{end}: start/end of input range relative
          to the "query space";
    \item \code{from_hit}: index of input range involved in hit;
    \item \code{to_hit}: index of (cigar,lmmpos) pair involved in hit.
  }
  The 4 columns are integer vectors.
}

\author{Valerie Obenchain and Hervé Pagès}

\seealso{
  \itemize{
    \item The \code{\link[GenomicAlignments]{mapToAlignments}} methods
          defined in the \pkg{GenomicAlignments} package.

    \item \code{\link{ref_pos_as_query_pos}} to project positions that
          are defined along the "reference space" onto the "query space".

    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \code{\link[IRanges]{findOverlaps}()} in the \pkg{IRanges} package.
  }
}

\examples{
set.seed(888)

## Random input ranges:
start <- sample(50000L, 10000, replace=TRUE)
end <- start + sample(15L, 10000, replace=TRUE) - 1L

## Random (cigar,lmmpos) pairs, kind of:
cigars <- sample(c("4M", "5M3I4M", "4M3D5M", "3M", "10M", "5M8N5M"),
                 25000, replace=TRUE)
lmmpos <- sample(50000L, 25000, replace=TRUE)

## map_ref_ranges_to_query():
system.time(df <- map_ref_ranges_to_query(start, end, cigars, lmmpos))
dim(df)
df[1:15, ]

## fast_map_ref_ranges_to_query() is about 300x-400x faster:
system.time(df2 <- fast_map_ref_ranges_to_query(start, end, cigars, lmmpos,
                                                strictly.sort.hits=TRUE))
stopifnot(identical(df, df2))
}

\keyword{manip}
