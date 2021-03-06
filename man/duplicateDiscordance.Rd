\name{duplicateDiscordance}
\alias{duplicateDiscordance}
\alias{duplicateDiscordance,SeqVarData,missing-method}
\alias{duplicateDiscordance,SeqVarData,SeqVarData-method}
\alias{duplicateDiscordance,SeqVarIterator,missing-method}

\title{Duplicate discordance}

\description{Find discordance rate for duplicate sample pairs}

\usage{
%\S4method{duplicateDiscordance}{SeqVarData,missing}(gdsobj, match.samples.on="subject.id", check.phase=FALSE, verbose=TRUE)
\S4method{duplicateDiscordance}{SeqVarData,missing}(gdsobj, match.samples.on="subject.id", by.variant=FALSE,
    all.pairs=TRUE, verbose=TRUE)
\S4method{duplicateDiscordance}{SeqVarIterator,missing}(gdsobj, match.samples.on="subject.id", by.variant=FALSE,
    all.pairs=TRUE, verbose=TRUE)
\S4method{duplicateDiscordance}{SeqVarData,SeqVarData}(gdsobj, obj2, match.samples.on=c("subject.id", "subject.id"),
    match.variants.on=c("alleles", "position"),
    discordance.type=c("genotype", "hethom"),
    by.variant=FALSE, verbose=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarData}} object with VCF data.}
  \item{obj2}{A \code{\link{SeqVarData}} object with VCF data.}
  \item{match.samples.on}{Character string or vector of strings indicating which column should be used for matching samples. See details.}
  \item{match.variants.on}{Character string of length one indicating how to match variants. See details.}
  \item{discordance.type}{Character string describing how discordances should be calculated. See details.}
  %\item{check.phase}{A logical indicating whether phase should be considered when calculating discordance.}
  \item{by.variant}{Calculate discordance by variant, otherwise by sample}
  \item{all.pairs}{Logical for whether to include all possible pairs of samples (\code{all.pairs=TRUE}) or only the first pair per subject (\code{all.pairs=FALSE}).}
  \item{verbose}{A logical indicating whether to print progress messages.}
}

\details{
  %For calls that involve only one gds file, duplicate discordance is calculated by sample pair and by variant.  If there are more than two samples per subject in \code{samples}, only the first two samples are used and a warning message is printed.
  %If \code{check.phase=TRUE}, variants with mismatched phase are considered discordant.  If \code{check.phase=FALSE}, phase is ignored.
  For calls that involve only one gds file, duplicate discordance is calculated by matching samples on common values of a column in \code{\link{sampleData}}. If \code{all.pairs=TRUE}, every possible pair of samples is included, so there may be multiple pairs per subject. If \code{all.pairs=FALSE}, only the first pair for each subject is used.
 
  For calls that involve two gds files, 
  duplicate discordance is calculated by matching sample pairs and variants between the two data sets. 
  Only biallelic SNVs are considered in the comparison.
  Variants can be matched using chromosome and position only (\code{match.variants.on="position"}) or by using chromosome, position, and alleles (\code{match.variants.on="alleles"}).
  If matching on alleles and the reference allele in the first dataset is the alternate allele in the second dataset, the genotype dosage will be recoded so the same allele is counted before making the comparison.
  If a variant in one dataset maps to multiple variants in the other dataset, only the first pair is considered for the comparison.
  Discordances can be calculated using either genotypes (\code{discordance.type = "genotype"}) or heterozygote/homozygote status (\code{discordance.type = "hethom"}).
  The latter is a method to calculate discordance that does not require alleles to be measured on the same strand in both datasets, so it is probably best to also set \code{match.variants.on = "position"} if using the \code{"hethom"} option.

  The argument \code{match.samples.on} can be used to select which column in the \code{\link{sampleData}} of the input \code{\link{SeqVarData}} object should be used for matching samples.
  For one gds file, \code{match.samples.on} should be a single string.
  For two gds files, \code{match.samples.on} should be a length-2 vector of character strings, where the first element is the column to use for the first gds object and the second element is the column to use for the second gds file.

  To exclude certain variants or samples from the calculate, use \code{\link{seqSetFilter}} to set appropriate filters on each gds object.

}

\value{
  %% For calls involving one gds file, a list with the following elements:
  %% \item{by.variant}{A data.frame with the number of discordances for
  %% each variant, the
  %% number of sample pairs with non-missing data, and the discordance rate
  %% (num.discord / num.pair).  Row names are variant ids.}
  %% \item{by.subject}{A data.frame with the sample ids for each pair, the
  %% number of discordances, the number of non-missing variants, and the
  %% discordance rate (num.discord / num.var).  Row.names are subject.id
  %% (as given in \code{samples}).}
  
  %% For calls involving two gds files,
  A data frame with the following columns, depending on whether \code{by.variant=TRUE} or \code{FALSE}:
\item{subject.id}{currently, this is the sample ID (\code{by.variant=FALSE} only)}
\item{sample.id.1/variant.id.1}{sample id or variant id in the first gds file}
\item{sample.id.2/variant.id.2}{sample id or variant id in the second gds file}
\item{n.variants/n.samples}{the number of non-missing variants or samples that were compared}
\item{n.concordant}{the number of concordant variants}
\item{n.alt}{the number of variants involving the alternate allele in either sample}
\item{n.alt.conc}{the number of concordant variants invovling the alternate allele in either sample}
\item{n.het.ref}{the number of mismatches where one call is a heterozygote and the other is a reference homozygote}
\item{n.het.alt}{the number of mismatches where one call is a heterozygote and the other is an alternate homozygote}
\item{n.ref.alt}{the number of mismatches where the calls are opposite homozygotes}

}

\author{Stephanie Gogarten, Adrienne Stilp}

\seealso{
  \code{\link{SeqVarData}},
  \code{\link{SeqVarIterator}}
}

\examples{
require(Biobase)

gds <- seqOpen(seqExampleFileName("gds"))

## the example file has one sample per subject, but we
## will match the first four samples into pairs as an example
sample.id <- seqGetData(gds, "sample.id")
samples <- AnnotatedDataFrame(data.frame(data.frame(subject.id=rep(c("subj1", "subj2"), times=45),
                      sample.id=sample.id,
                      stringsAsFactors=FALSE)))
seqData <- SeqVarData(gds, sampleData=samples)

## set a filter on the first four samples
seqSetFilter(seqData, sample.id=sample.id[1:4])

disc <- duplicateDiscordance(seqData, by.variant=FALSE)
disc
disc <- duplicateDiscordance(seqData, by.variant=TRUE)
head(disc)

## recommended to use an iterator object for large datasets
iterator <- SeqVarBlockIterator(seqData)
disc <- duplicateDiscordance(iterator, by.variant=FALSE)
disc

seqClose(gds)
}

\keyword{manip}
