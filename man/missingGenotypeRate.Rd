\name{missingGenotypeRate}
\alias{missingGenotypeRate}
\alias{missingGenotypeRate,SeqVarGDSClass-method}

\title{Missing genotype rate}

\description{
  Calculate missing genotype rate by variant or by sample
}

\usage{
\S4method{missingGenotypeRate}{SeqVarGDSClass}(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{margin}{Possible values are "by.variant" or "by.sample,"
    indicating whether the calculation should be done over all samples
    for each variant, or over all variants for each sample.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
}

\details{
  Calculates the fraction of missing genotypes in a GDS object, either
  by variant or by sample.
}

\value{
  A numeric vector of missing genotype rates.  If
  \code{margin="by.variant"}, the vector will have length equal to the
  number of variants in the GDS object.  If \code{margin="by.sample"},
  the vector will have length equal to the number of samples.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{getGenotype}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
head(missingGenotypeRate(gds, margin="by.variant"))
head(missingGenotypeRate(gds, margin="by.sample"))
seqClose(gds)
}

\keyword{manip}
