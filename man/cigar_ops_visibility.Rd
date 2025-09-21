\name{cigar_ops_visibility}

\alias{CIGAR_OPS}
\alias{cigar_ops_visibility}

\title{Visibility of CIGAR operations}

\description{
  CIGAR operations and their visibility in various \emph{projection spaces}.
}

\usage{
CIGAR_OPS

cigar_ops_visibility(ops=CIGAR_OPS)
}

\arguments{
  \item{ops}{
    Character vector where the elements are single letters representing
    valid CIGAR operations. Must be a subset of \code{CIGAR_OPS}.
  }
}

\details{
  The 8 supported \emph{projection spaces} are: "reference",
  "reference-N-regions-removed", "query", "query-before-hard-clipping",
  "query-after-soft-clipping", "pairwise", "pairwise-N-regions-removed",
  and "pairwise-dense".

  Each space can be characterized by the extended CIGAR
  operations that are \emph{visible} in it. A CIGAR operation is
  said to be \emph{visible} in a given space if it "runs along it",
  that is, if it's associated with a block of contiguous positions in
  that space (the size of the block being the length of the operation).
  For example, the M/=/X operations are \emph{visible} in all spaces,
  the D/N operations are \emph{visible} in the "reference" space
  but not in the "query" space, the S operation is \emph{visible}
  in the "query" space but not in the "reference" or in
  the "query-after-soft-clipping" space, etc...

  Here are the extended CIGAR operations that are \emph{visible}
  in each space:
  \enumerate{
    \item reference: M, D, N, =, X

    \item reference-N-regions-removed: M, D, =, X

    \item query: M, I, S, =, X

    \item query-before-hard-clipping: M, I, S, H, =, X

    \item query-after-soft-clipping: M, I, =, X

    \item pairwise: M, I, D, N, =, X

    \item pairwise-N-regions-removed: M, I, D, =, X

    \item pairwise-dense: M, =, X
  }
  Note that CIGAR operations M, =, and X are visible in all spaces.
}

\value{
  \code{CIGAR_OPS} is a predefined character vector containing the valid
  (extended) CIGAR operations: M, I, D, N, S, H, P, =, X.
  See official SAM/BAM Format specs at
  \url{https://samtools.github.io/hts-specs/SAMv1.pdf} for the list
  of extended CIGAR operations and their meaning.

  \code{cigar_ops_visibility()} returns an 8-row integer matrix with 1 row
  per space and 1 column per CIGAR operation. The matrix is made of 0's
  and 1's indicating visibility.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{explode_cigars} to extract the letters (or lengths) of
          the CIGAR operations contained in a vector of CIGAR strings.

    \item \code{\link{tabulate_cigar_ops}} to count the occurences of CIGAR
           operations in a vector of CIGAR strings.

    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \code{\link{trim_cigars_along_ref}} and
          \code{\link{trim_cigars_along_query}} to trim CIGAR strings
          along the "reference space" and "query space", respectively.

    \item \link{cigars_as_ranges} to turn CIGAR strings into ranges
          of positions.

    \item \code{\link{project_positions}} to project positions from query
          to reference space and vice versa.

    \item \code{\link{project_sequences}} to project sequences from one
          space to the other.
  }
}

\examples{
CIGAR_OPS  # valid CIGAR operations

cigar_ops_visibility()  # visibility in each "projection space"
}

\keyword{manip}
