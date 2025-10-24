\name{cigarillo-package}

\alias{cigarillo-package}
\alias{cigarillo}
\docType{package}

\title{Efficient manipulation of CIGAR strings}

\description{
  CIGAR stands for Concise Idiosyncratic Gapped Alignment Report.
  CIGAR strings are found in the BAM files produced by most aligners
  and in the AIRR-formatted output produced by IgBLAST.

  The \pkg{cigarillo} package provides functions to parse and inspect CIGAR
  strings, trim them, turn them into ranges of positions relative to
  the "query space" or "reference space", and project positions or
  sequences from one space to the other. Note that these operations
  are low-level operations that the user rarely needs to perform
  directly. More typically, they are performed behind the scene by
  higher-level functionality implemented in other packages like
  Bioconductor packages \pkg{GenomicAlignments} and \pkg{igblastr}.
}

\details{
  For an overview of the functionality provided by the package, please
  refer to the vignette:
  \preformatted{    vignette("cigarillo", package="cigarillo")}
}

\author{
  Hervé Pagès, Valerie Obenchain, Michael Lawrence

  With contributions from Martin Morgan, Patrick Aboyoun, and Fedor Bezrukov
}

\seealso{
  \itemize{
    \item \code{\link{cigar_ops_visibility}} for an introduction to CIGAR
          operations and their visibility in various "projection spaces".

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

    \item The \link[IRanges]{RleList} class in the \pkg{IRanges} package.
  }
}

\keyword{package}
