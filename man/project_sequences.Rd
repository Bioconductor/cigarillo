\name{project_sequences}

\alias{project_sequences}


\title{Project sequences from one space to the other}

\description{
  \code{project_sequences} can project sequences that belong to a given
  space (e.g. the "query space") onto another space (e.g. the
  "reference space") by removing/injecting substrings from/into them,
  based on the associated CIGAR strings.

  Its primary use case is to project the read sequences stored in a
  BAM file (which are considered to belong to the "query space") onto
  the "reference space". It can also be used to remove the parts of the
  read sequences that correspond to soft-clipping. More generally it can
  project sequences that belong to any supported space onto any other
  supported space. See the Details section below for the list of supported
  spaces.
}

\usage{
project_sequences(x, cigars, from="query", to="reference",
                  I.letter="-", D.letter="-", N.letter=".",
                  S.letter="+", H.letter="+")
}

\arguments{
  \item{x}{
    An \link[Biostrings]{XStringSet} derivative (e.g.
    \link[Biostrings]{BStringSet}, \link[Biostrings]{DNAStringSet},
    or \link[Biostrings]{AAStringSet} object) containing sequences
    that are considered to belong to the \code{from} space (see below).
  }
  \item{cigars}{
    A character vector (or factor) parallel to \code{x}
    containing CIGAR strings.
  }
  \item{from, to}{
    A single string specifying one of the 8 supported "projection spaces".
    See \code{?\link{cigar_ops_visibility}} for more information.
    \code{from} must be the current space (i.e. the space that
    the sequences in \code{x} belong to) and \code{to} is the
    space onto which the sequences in \code{x} must be projected.
  }
  \item{I.letter, D.letter, N.letter, S.letter, H.letter}{
    A single letter used as a filler for injections. More on this in
    the Details section below.
  }
}

\details{
  \code{project_sequences} projects a sequence that belongs to one space
  onto another by (1) removing the substrings associated with operations
  that are no longer \emph{visible} in the new space, and (2) injecting
  substrings associated with operations that become \emph{visible} in the
  new space. Each injected substring has the length of the operation
  associated with it, and its content is controlled via the corresponding
  \code{*.letter} argument.

  For example, when going from the "query" space to the "reference" space
  (the default), the I- and S-substrings (i.e. the substrings associated
  with I/S operations) are removed, and substrings associated with D/N
  operations are injected. More precisely, the D-substrings are filled
  with the letter specified in \code{D.letter}, and the N-substrings with
  the letter specified in \code{N.letter}. The other \code{*.letter}
  arguments are ignored in that case.
}

\value{
  An \link{XStringSet} derivative of the same class as input object
  \code{x}, and parallel to \code{x}. The names on \code{x}, if any,
  are propagated.
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

    \item \link{cigars_as_ranges} to turn CIGAR strings into ranges
          of positions.

    \item \code{\link{project_positions}} to project positions from query
          to reference space and vice versa.

    \item The \code{\link[GenomicAlignments]{stackStringsFromBam}} function
          in the \pkg{GenomicAlignments} package for stacking the read
          sequences (or their quality strings) stored in a BAM file on a
          region of interest.

    \item The \code{\link[GenomicAlignments]{readGAlignments}} function
          in the \pkg{GenomicAlignments} package for loading read sequences
          from a BAM file (as a \link[GenomicAlignments]{GAlignments} object).

    \item The \code{\link[Biostrings]{extractAt}} and
          \code{\link[Biostrings]{replaceAt}} functions in the \pkg{Biostrings}
          package for extracting/replacing arbitrary substrings from/in a
          string or set of strings.
  }
}

\examples{
library(GenomicAlignments)

## ---------------------------------------------------------------------
## A. FROM "query" TO "reference" SPACE
## ---------------------------------------------------------------------

## Load read sequences from a BAM file (they will be returned in a
## GAlignments object):
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what="seq")
gal <- readGAlignments(bamfile, param=param)
qseq <- mcols(gal)$seq  # the read sequences (aka query sequences)

## Project the query sequences onto the reference space. This will
## remove the substrings associated with insertions to the reference
## (I operations) and soft clipping (S operations), and will inject new
## substrings (filled with "-") where deletions from the reference (D
## operations) and skipped regions from the reference (N operations)
## occurred during the alignment process:
qseq_on_ref <- project_sequences(qseq, cigar(gal))

## A typical use case for doing the above is to compute 1 consensus
## sequence per chromosome. The code below shows how this can be done
## in 2 extra steps.

## Step 1: Compute one consensus matrix per chromosome.
qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
pos_by_chrom <- splitAsList(start(gal), seqnames(gal))

cm_by_chrom <- lapply(names(pos_by_chrom),
    function(seqname)
        consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                        as.prob=TRUE,
                        shift=pos_by_chrom[[seqname]]-1,
                        width=seqlengths(gal)[[seqname]]))
names(cm_by_chrom) <- names(pos_by_chrom)

## 'cm_by_chrom' is a list of consensus matrices. Each matrix has 17
## rows (1 per letter in the DNA alphabet) and 1 column per chromosome
## position.

## Step 2: Compute the consensus string from each consensus matrix.
## We'll put "+" in the strings wherever there is no coverage for that
## position, and "N" where there is coverage but no consensus.
cs_by_chrom <- lapply(cm_by_chrom,
    function(cm) {
        ## Because consensusString() doesn't like consensus matrices
        ## with columns that contain only zeroes (and you will have
        ## columns like that for chromosome positions that don't
        ## receive any coverage), we need to "fix" 'cm' first.
        idx <- colSums(cm) == 0
        cm["+", idx] <- 1
        DNAString(consensusString(cm, ambiguityMap="N"))
    })

## consensusString() provides some flexibility to let you extract
## the consensus in different ways. See '?consensusString' in the
## Biostrings package for the details.

## Finally, note that the read quality strings can also be used as
## input for project_sequences():
param <- ScanBamParam(what="qual")
gal <- readGAlignments(bamfile, param=param)
qual <- mcols(gal)$qual  # the read quality strings
qual_on_ref <- project_sequences(qual, cigar(gal))
## Note that since the "-" letter is a valid quality code, there is
## no way to distinguish it from the "-" letters inserted by
## project_sequences().

## ---------------------------------------------------------------------
## B. FROM "query" TO "query-after-soft-clipping" SPACE
## ---------------------------------------------------------------------

## Going from "query" to "query-after-soft-clipping" simply removes
## the substrings associated with soft clipping (S operations):
qseq <- DNAStringSet(c("AAAGTTCGAA", "TTACGATTAN", "GGATAATTTT"))
cigars <- c("3H10M", "2S7M1S2H", "2M1I1M3D2M4S")
clipped_qseq <- project_sequences(qseq, cigars,
                                  from="query",
                                  to="query-after-soft-clipping")

project_sequences(clipped_qseq, cigars,
                  from="query-after-soft-clipping", to="query")

project_sequences(clipped_qseq, cigars,
                  from="query-after-soft-clipping", to="query",
                  S.letter="-")

## ---------------------------------------------------------------------
## C. BRING QUERY AND REFERENCE SEQUENCES TO THE "pairwise"
##    OR "pairwise-dense" SPACE
## ---------------------------------------------------------------------

## Load read sequences from a BAM file:
library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
param <- ScanBamParam(what="seq",
                      which=GRanges("chr14", IRanges(1, 25000000)))
gal <- readGAlignments(bamfile, param=param)
qseq <- mcols(gal)$seq  # the read sequences (aka query sequences)

## Load the corresponding reference sequences from the appropriate
## BSgenome package (the reads in RNAseqData.HNRNPC.bam.chr14 were
## aligned to hg19):
library(BSgenome.Hsapiens.UCSC.hg19)
rseq <- getSeq(Hsapiens, as(gal, "GRanges"))  # the reference sequences

## Bring 'qseq' and 'rseq' to the "pairwise" space.
## For 'qseq', this will remove the substrings associated with soft
## clipping (S operations) and inject substrings (filled with "-")
## associated with deletions from the reference (D operations) and
## skipped regions from the reference (N operations). For 'rseq', this
## will inject substrings (filled with "-") associated with insertions
## to the reference (I operations).
qseq2 <- project_sequences(qseq, cigar(gal),
                           from="query", to="pairwise")
rseq2 <- project_sequences(rseq, cigar(gal),
                           from="reference", to="pairwise")

## Sanity check: 'qseq2' and 'rseq2' should have the same shape.
stopifnot(identical(elementNROWS(qseq2), elementNROWS(rseq2)))

## A closer look at reads with insertions and deletions:
cigar_op_table <- cigarOpTable(cigar(gal))
head(cigar_op_table)

I_idx <- which(cigar_op_table[ , "I"] >= 2)  # at least 2 insertions
qseq2[I_idx]
rseq2[I_idx]

D_idx <- which(cigar_op_table[ , "D"] >= 2)  # at least 2 deletions
qseq2[D_idx]
rseq2[D_idx]

## A closer look at reads with skipped regions:
N_idx <- which(cigar_op_table[ , "N"] != 0)
qseq2[N_idx]
rseq2[N_idx]

## A variant of the "pairwise" space is the "pairwise-dense" space.
## In that space, all indels and skipped regions are removed from 'qseq'
## and 'rseq'.
qseq3 <- project_sequences(qseq, cigar(gal),
                           from="query", to="pairwise-dense")
rseq3 <- project_sequences(rseq, cigar(gal),
                           from="reference", to="pairwise-dense")

## Sanity check: 'qseq3' and 'rseq3' should have the same shape.
stopifnot(identical(elementNROWS(qseq3), elementNROWS(rseq3)))

## Insertions were removed:
qseq3[I_idx]
rseq3[I_idx]

## Deletions were removed:
qseq3[D_idx]
rseq3[D_idx]

## Skipped regions were removed:
qseq3[N_idx]
rseq3[N_idx]
}

\keyword{methods}
\keyword{manip}
