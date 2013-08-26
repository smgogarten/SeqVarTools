\name{duplicateDiscordance}
\alias{duplicateDiscordance}
\alias{duplicateDiscordance,SeqVarGDSClass-method}

\title{Duplicate discordance}

\description{Find discordance rate for duplicate sample pairs}

\usage{
\S4method{duplicateDiscordance}{SeqVarGDSClass}(gdsobj, samples, check.phase=FALSE, verbose=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{samples}{A data.frame with columns (sample.id, subject.id).
    "sample.id" values should correspond to "sample.id" in \code{gdsobj}.
    "subject.id" should match for duplicate samples.}
  \item{check.phase}{A logical indicating whether phase should be
  considered when calculating discordance.}
  \item{verbose}{A logical indicating whether to print a progress
  message for each sample.}
}

\details{
  Duplicate discordance is calculated by sample pair and by variant.  If
  there are more than two samples per subject in \code{samples}, only the first two
  samples are used and a warning message is printed.

  If \code{check.phase=TRUE}, variants with mismatched phase are
  considered discordant.  If \code{check.phase=FALSE}, phase is ignored.
}

\value{
  A list with the following elements:
  \item{by.variant}{A data.frame with the number of discordances for
  each variant, the
  number of sample pairs with non-missing data, and the discordance rate
  (num.discord / num.pair).  Row names are variant ids.}
  \item{by.subject}{A data.frame with the sample ids for each pair, the
  number of discordances, the number of non-missing variants, and the
  discordance rate (num.discord / num.var).  Row.names are subject.id
  (as given in \code{samples}).}
}

\author{Stephanie Gogarten}

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
disc <- duplicateDiscordance(gds, samples)
head(disc$by.variant)
disc$by.subject
seqClose(gds)
}

\keyword{manip}
