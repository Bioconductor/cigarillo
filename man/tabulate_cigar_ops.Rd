\name{tabulate_cigar_ops}

\alias{tabulate_cigar_ops}

\title{Tabulate CIGAR operations}

\description{
  Count the occurences of CIGAR operations in a vector of CIGAR strings.
}

\usage{
tabulate_cigar_ops(cigars, oplens.as.weights=FALSE)
}

\arguments{
  \item{cigars}{
    A character vector (or factor) containing CIGAR strings.
  }
  \item{oplens.as.weights}{
    \code{TRUE} or \code{FALSE}.

    Should the operation lengths be used as weights for the counts?
  }
}

\value{
  An integer matrix with 1 row per CIGAR string in \code{cigars}
  and 1 column per CIGAR operation in \code{CIGAR_OPS}.
}

\author{Patrick Aboyoun and Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{cigar_ops_visibility}} for an introduction to CIGAR
          operations and their visibility in various "projection spaces".

    \item \link{explode_cigars} to extract the letters (or lengths) of
          the CIGAR operations contained in a vector of CIGAR strings.

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
my_cigars <- c(
    "40M2I9M",
    "60M",
    "3H15M55N4M2I6M2D5M6S",
    "50=2X3=1X10=",
    "2S10M2000N15M",
    "3H33M5H"
)

op_counts <- tabulate_cigar_ops(my_cigars)
op_counts

tabulate_cigar_ops(my_cigars, oplens.as.weights=TRUE)

## Get the total number of operations per CIGAR string:
rowSums(op_counts)  # a numeric vector parallel to 'my_cigars'

## Note that the above is equivalent to -- but much faster and more
## memory-efficient than -- 'lengths(explode_cigar_ops(my_cigars))'
## or 'lengths(explode_cigar_oplens(my_cigars))':
nop_per_cig <- as.integer(rowSums(op_counts))
stopifnot(
    identical(nop_per_cig, lengths(explode_cigar_ops(my_cigars))),
    identical(nop_per_cig, lengths(explode_cigar_oplens(my_cigars)))
)

## Identify CIGAR strings with indels:
has_indels <- rowSums(op_counts[ , c("I", "D")]) != 0
has_indels  # a logical vector parallel to 'my_cigars'

## Summarize the counts for the whole vector of CIGAR strings:
colSums(op_counts)
}

\keyword{manip}
