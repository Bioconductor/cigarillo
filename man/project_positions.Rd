\name{project_positions}

\alias{project_positions}

\alias{query_pos_as_ref_pos}
\alias{ref_pos_as_query_pos}

\title{Project positions from query to reference space and vice versa}

\description{
  \code{query_pos_as_ref_pos()} projects positions defined along
  the "query space" onto the "reference space", that is, it turns
  them into positions defined along the "reference space".

  \code{ref_pos_as_query_pos()} does the opposite i.e. it projects
  positions that are defined along the "reference space" onto
  the "query space".
}

\usage{
query_pos_as_ref_pos(query_pos, cigars, lmmpos, narrow.left)

ref_pos_as_query_pos(ref_pos, cigars, lmmpos, narrow.left)
}

\arguments{
  \item{query_pos}{
    An integer vector containing positions relative to the "query space".
  }
  \item{cigars}{
    A character vector (or factor) parallel to \code{query_pos}
    containing CIGAR strings.
  }
  \item{lmmpos}{
    An integer vector parallel to \code{cigars} and \code{query_pos}.
    For each CIGAR string in \code{cigars}, \code{lmmpos} must contain
    the 1-based leftmost mapping POSition of the alignment described
    by the CIGAR string. Note that these positions must be relative to
    the "reference space".
  }
  \item{ref_pos}{
    An integer vector containing positions relative to the "reference space".
  }
  \item{narrow.left}{
    For \code{query_pos_as_ref_pos()}: How should positions in
    the "query space" that fall within an insertion be treated?

    Such positions are peculiar, because, strictly speaking, they don't
    have corresponding positions in the "reference space". Instead,
    each of them falls \emph{between} two adjacent positions in
    the "reference space". Another way to describe this situation is
    to say that each of them is mapped to a zero-width range along
    the "reference space".

    If \code{narrow.left} is \code{TRUE}, such position will be mapped
    to the position that is immediately on the left of the corresponding
    zero-width range on the "reference space".
    If \code{narrow.left} is \code{FALSE}, it will be mapped to the
    position that is immediately on the right of the corresponding
    zero-width range on the "reference space".

    For \code{ref_pos_as_query_pos()}: How should positions in
    the "reference space" that fall within a deletion be treated?

    Such positions are peculiar, because, strictly speaking, they don't
    have corresponding positions in the "query space". Instead,
    each of them falls \emph{between} two adjacent positions in
    the "query space". Another way to describe this situation is
    to say that each of them is mapped to a zero-width range along
    the "query space".

    If \code{narrow.left} is \code{TRUE}, such position will be mapped
    to the position that is immediately on the left of the corresponding
    zero-width range on the "query space".
    If \code{narrow.left} is \code{FALSE}, it will be mapped to the
    position that is immediately on the right of the corresponding
    zero-width range on the "query space".
  }
}

\value{
  An integer vector parallel to the input positions. \code{NA}s in
  the returned vector indicate input positions that cannot be mapped.
}

\author{Michael Lawrence}

\seealso{
  \itemize{
    \item \link{inspect_cigars} for basic utility functions to parse
          and inspect CIGAR strings.

    \item \link{cigar_extent} for functions that calculate the \emph{extent}
          of a CIGAR string, that is, the number of positions spanned by
          the alignment that it describes.

    \item \code{\link{trim_cigars_along_ref}} and
          \code{\link{trim_cigars_along_query}} to trim CIGAR strings
          along the "reference space" and "query space", respectively.

    \item \link{cigars_as_ranges} to turn CIGAR strings into ranges
          of positions.

    \item \link{project_sequences} to project sequences from one space
          to the other.
  }
}

\examples{
query_pos <- -1:11
cigars <- rep("5M3I2M", 13)
lmmpos <- rep(101, 13)
query_pos_as_ref_pos(query_pos, cigars, lmmpos, narrow.left=TRUE)
query_pos_as_ref_pos(query_pos, cigars, lmmpos, narrow.left=FALSE)
}

\keyword{manip}
