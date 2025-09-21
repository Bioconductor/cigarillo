\name{explode_cigars}

\alias{explode_cigars}

\alias{validate_cigars}

\alias{explode_cigar_ops}
\alias{explode_cigar_oplens}
\alias{cigars_as_RleList}

\title{Explode CIGAR strings}

\description{
  Use \code{explode_cigar_ops()} (or \code{explode_cigar_oplens()}) to
  extract the letters (or lengths) of the CIGAR operations contained
  in a vector of CIGAR strings.
}

\usage{
explode_cigar_ops(cigars, ops=CIGAR_OPS)
explode_cigar_oplens(cigars, ops=CIGAR_OPS)

cigars_as_RleList(cigars)
}

\arguments{
  \item{cigars}{
    A character vector (or factor) containing CIGAR strings.
  }
  \item{ops}{
    Character vector where the elements are single letters representing
    valid CIGAR operations. Must be a subset of \code{\link{CIGAR_OPS}}.
    See \code{?\link{CIGAR_OPS}} for more information.

    \code{explode_cigar_ops()} and \code{explode_cigar_oplens()} will
    ignore operations not listed in \code{ops} (in addition to 0-length
    operations which are always ignored).
  }
}

\value{
  For \code{explode_cigar_ops} and \code{explode_cigar_oplens}:
  Both functions return a list parallel to \code{cigars} where each
  list element is a character vector (for \code{explode_cigar_ops})
  or an integer vector (for \code{explode_cigar_oplens}). The two lists
  are guaranteed to have the same shape, that is, the same \code{length()}
  and same \code{lengths()}.

  More precisely: The i-th character vector in the list returned by
  \code{explode_cigar_ops} contains one single-letter string per CIGAR
  operation in \code{cigars[i]}. The i-th integer vector in the list
  returned by \code{explode_cigar_oplens} contains the corresponding
  CIGAR operation lengths. Operations not listed in \code{ops} and 0-length
  operations are ignored.

  For \code{cigars_as_RleList}: An \link[IRanges]{RleList} object.
}

\author{Hervé Pagès, Martin Morgan, and Patrick Aboyoun}

\seealso{
  \itemize{
    \item \code{\link{cigar_ops_visibility}} for an introduction to CIGAR
          operations and their visibility in various "projection spaces".

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

    \item The \link[IRanges]{RleList} class in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## Turn CIGAR strings into other useful representations
## ---------------------------------------------------------------------

my_cigars <- c(
    "40M2I9M",
    "60M",
    "3H15M55N4M2I6M2D5M6S",
    "50=2X3=1X10=",
    "2S10M2000N15M",
    "3H33M5H"
)

cig_ops <- explode_cigar_ops(my_cigars)
cig_ops

cig_oplens <- explode_cigar_oplens(my_cigars)
cig_oplens

explode_cigar_ops(my_cigars, ops=c("I", "S"))
explode_cigar_oplens(my_cigars, ops=c("I", "S"))

cigs_as_rlelist <- cigars_as_RleList(my_cigars)
cigs_as_rlelist

## ---------------------------------------------------------------------
## Results can be coerced to CharacterList or IntegerList
## ---------------------------------------------------------------------

as(cig_ops, "CharacterList")
as(cig_oplens, "IntegerList")
as(cigs_as_rlelist, "CharacterList")

## ---------------------------------------------------------------------
## Sanity checks
## ---------------------------------------------------------------------

stopifnot(
    identical(as.list(runValue(cigs_as_rlelist)), cig_ops),
    identical(as.list(runLength(cigs_as_rlelist)), cig_oplens)
)
}

\keyword{manip}
