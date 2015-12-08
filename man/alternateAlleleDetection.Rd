\name{alternateAlleleDetection}
\alias{alternateAlleleDetection}
\alias{alternateAlleleDetection,SeqVarData,SeqVarData-method}

\title{alternateAlleleDetection}

\description{Calculate rates of detecting minor alleles given a ``gold standard'' dataset}

\usage{
\S4method{alternateAlleleDetection}{SeqVarData,SeqVarData}(gdsobj, gdsobj2, match.samples.on=c("subject.id", "subject.id"), verbose=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarData}} object with VCF data.}
  \item{gdsobj2}{A \code{\link{SeqVarData}} object with VCF data to be used as the ``gold standard''.}
  \item{match.samples.on}{A length-2 character vector indicating the column to be used for matching in each dataset's \code{sampleData} annotation}
  \item{verbose}{A logical indicating whether to print progress messages.}
}

\details{

  Calculates the accuracy of detecting alternate alleles in one dataset (\code{gdsobj}) given a ``gold standard'' dataset (\code{gdsobj2}).
  Samples are matched using the \code{match.samples.on} argument. The first element of \code{match.samples.on} indicates the column to be used as the subject identifier for the first dataset, and the second element is the column to be used for the second dataset.
  Variants are matched on position and alleles using bi-allelic SNVs only.
  No genotype recoding is done; the variant must have the same reference and alternate alleles to be considered a match.
  If a variant in one dataset matches to multiple variants in the second dataset, then only the first match will be used.
  If a variant is missing in either dataset for a given sample pair, that sample pair is ignored for that variant.
  To exclude certain variants or samples from the calculate, use \code{\link{seqSetFilter}} to set appropriate filters on each gds object.
  
  This test is positive if an alternate allele was been detected.
  Results are returned on an allele level, such that:

  \code{TP}, \code{TN}, \code{FP}, and \code{FN} are calculated as follows:
  \tabular{lllll}{
              \tab    \tab           \tab genoData2 \tab \cr
              \tab    \tab aa        \tab Ra        \tab RR\cr
              \tab aa \tab 2TP       \tab 1TP + 1FP \tab 2FP\cr
    genoData1 \tab Ra \tab 1TP + 1FN \tab 1TN + 1TP \tab 1TN + 1FP\cr
	      \tab RR \tab 2FN       \tab 1FN + 1TN \tab 2TN\cr
	    }
  where ``R'' indicates a reference allele and ``a'' indicates an alternate allele.
}

\value{
  A data frame with the following columns:
\item{variant.id.1}{variant id from the first dataset}
\item{variant.id.2}{matched variant id from the second dataset}
\item{n.samples}{the number of samples with non-missing data for this variant}
\item{true.pos}{the number of alleles that are true positives for this variant}
\item{true.neg}{the number of alleles that are true negatives for this variant}
\item{false.pos}{the number of alleles that are false positives for this variant}
\item{false.neg}{the number of alleles that are false negatives for this variant}
}

\author{Adrienne Stilp}

\seealso{
  \code{\link{SeqVarGDSClass}}
}

\examples{
\dontrun{
gds1 <- seqOpen(gdsfile.1) # dataset to test, e.g. sequencing
sample1 <- data.frame(subject.id=c("a", "b", "c"), sample.id=c("A", "B", "C"), stringsAsFactors=F)
seqData1 <- SeqVarData(gds1, sampleData=sample1)

gds2 <- seqOpen(gdsfile.2) # gold standard dataset, e.g. array genotyping
sample2 <- data.frame(subject.id=c("b", "c", "d"), sample.id=c("B", "C", "D"), stringsAsFactors=F)
seqData2 <- SeqVarData(gds2, sampleData=sample2)

res <- alleleDetectionAccuracy(seqData1, seqData2)
}
}

\keyword{manip}
