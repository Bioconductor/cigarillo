\name{inspect-cigars}

\alias{inspect-cigars}
\alias{inspect_cigars}

\alias{validate_cigars}

\alias{CIGAR_OPS}
\alias{explode_cigar_ops}
\alias{explode_cigar_oplens}
\alias{cigars_as_RleList}
\alias{tabulate_cigar_ops}

\title{Inspect CIGAR strings}

\description{
  Basic utility functions to parse and inspect CIGAR strings.
}

\usage{
CIGAR_OPS

explode_cigar_ops(cigars, ops=CIGAR_OPS)
explode_cigar_oplens(cigars, ops=CIGAR_OPS)
cigars_as_RleList(cigars)

tabulate_cigar_ops(cigars)
}

\arguments{
  \item{cigars}{
    A character vector or factor containing CIGAR strings.
  }
  \item{ops}{
    Character vector containing the CIGAR operations to consider.
    This must be a subset of \code{CIGAR_OPS}. \code{explode_cigar_ops()}
    and \code{explode_cigar_oplens()} will ignore operations not listed
    in \code{ops} (in addition to 0-length operations which are always
    ignored).
  }
}

\value{
  \code{CIGAR_OPS} is a predefined character vector containing the valid
  (extended) CIGAR operations: M, I, D, N, S, H, P, =, X.
  See official SAM/BAM Format specs at
  \url{https://samtools.github.io/hts-specs/SAMv1.pdf} for the list of
  extended CIGAR operations and their meaning.

  For \code{explode_cigar_ops} and \code{explode_cigar_oplens}:
  Both functions return a list of the same length as \code{cigars} where
  each list element is a character vector (for \code{explode_cigar_ops})
  or an integer vector (for \code{explode_cigar_oplens}). The two lists
  are guaranteed to have the same shape, that is, the same \code{length()}
  and same \code{lengths()}.
  The i-th character vector in the list returned by \code{explode_cigar_ops}
  contains one single-letter string per CIGAR operation in \code{cigars[i]}.
  The i-th integer vector in the list returned by \code{explode_cigar_oplens}
  contains the corresponding CIGAR operation lengths. Operations not listed
  in \code{ops} and 0-length operations are ignored.

  For \code{cigars_as_RleList}: An \link[IRanges]{RleList} object.

  For \code{tabulate_cigar_ops}: An integer matrix with number of rows equal
  to the length of \code{cigars} and nine columns, one for each extended
  CIGAR operation.
}

\author{Hervé Pagès & P. Aboyoun}

\seealso{
  \itemize{
    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \code{\link{trim_cigars_along_ref}} and
          \code{\link{trim_cigars_along_query}} to trim CIGAR strings
          along the "reference space" and "query space", respectively.

    \item The \link[IRanges]{RleList} class in the \pkg{IRanges} package.
  }
}

\examples{
## Valid CIGAR operations:
CIGAR_OPS

## Turn CIGAR strings into other useful representations:
my_cigars <- c("40M2I9M", "3H15M55N4M2I6M2D5M6S",
               "2S10M2000N15M", "3H33M5H")

explode_cigar_ops(my_cigars)
explode_cigar_oplens(my_cigars)

explode_cigar_ops(my_cigars, ops=c("I", "S"))
explode_cigar_oplens(my_cigars, ops=c("I", "S"))

cigars_as_RleList(my_cigars)

## The RleList object returned by cigars_as_RleList() can be coerced
## to CharacterList:
as(cigars_as_RleList(my_cigars), "CharacterList")

## Tabulate CIGAR operations:
tabulate_cigar_ops(my_cigars)
}

\keyword{manip}
