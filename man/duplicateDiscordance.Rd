\name{duplicateDiscordance}
\alias{duplicateDiscordance}
\alias{duplicateDiscordance,SeqVarGDSClass,missing-method}
\alias{duplicateDiscordance,SeqVarGDSClass,SeqVarGDSClass-method}

\title{Duplicate discordance}

\description{Find discordance rate for duplicate sample pairs}

\usage{
\S4method{duplicateDiscordance}{SeqVarGDSClass,missing}(gdsobj, samples=NULL, check.phase=FALSE, verbose=TRUE)
\S4method{duplicateDiscordance}{SeqVarGDSClass,SeqVarGDSClass}(gdsobj, obj2, by.variant=FALSE, verbose=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{obj2}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{samples}{A data.frame with columns (sample.id, subject.id).
    "sample.id" values should correspond to "sample.id" in \code{gdsobj}.
    "subject.id" should match for duplicate samples.}
  \item{check.phase}{A logical indicating whether phase should be
  considered when calculating discordance.}
  \item{by.variant}{Calculate discordance by variant, otherwise by sample}
  \item{verbose}{A logical indicating whether to print progress messages.}
}

\details{
  For calls that involve only one gds file,
  duplicate discordance is calculated by sample pair and by variant.  If
  there are more than two samples per subject in \code{samples}, only the first two
  samples are used and a warning message is printed.
  If \code{check.phase=TRUE}, variants with mismatched phase are
  considered discordant.  If \code{check.phase=FALSE}, phase is ignored.
  
  For calls that invovle two gds files,
  duplicate discordance is calculated by matching sample pairs (not subjects).
  Samples are matched only if they have the same sample.id in each gds object, so the samples must have the same identifier in both gds files.
  In the future, support will be added for matching on a subjectID instead of on sample IDs.
  Variants are matched using position and alleles.
  Only biallelic SNVs are considered.
  No flipping is done; the same allele must be the reference in both gds files.
  
  To exclude certain variants or samples from the calculate, use \code{\link{seqSetFilter}} to set appropriate filters on each gds object.

}

\value{
  For calls involving one gds file, a list with the following elements:
  \item{by.variant}{A data.frame with the number of discordances for
  each variant, the
  number of sample pairs with non-missing data, and the discordance rate
  (num.discord / num.pair).  Row names are variant ids.}
  \item{by.subject}{A data.frame with the sample ids for each pair, the
  number of discordances, the number of non-missing variants, and the
  discordance rate (num.discord / num.var).  Row.names are subject.id
  (as given in \code{samples}).}
  
  For calls invovling two gds files,
  A data frame with the following columns, depending on whether \code{by.variant=TRUE} or \code{FALSE}:
\item{subjectID}{currently, this is the sample ID (\code{by.variant=FALSE} only)}
\item{sample.id.1/variant.id.1}{sample id or variant id in the first gds file}
\item{sample.id.2/variant.id.1}{sample id or variant id in the second gds file}
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
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
## the example file has one sample per subject, but we
## will match the first four samples into pairs as an example
sample.id <- seqGetData(gds, "sample.id")
samples <- data.frame(subject.id=rep(c("subj1", "subj2"), each=2),
                      sample.id=sample.id[1:4],
                      stringsAsFactors=FALSE)
disc <- duplicateDiscordance(gds, samples=samples)
head(disc$by.variant)
disc$by.subject
seqClose(gds)
}

\keyword{manip}
