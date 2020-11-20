\name{heterozygosity}
\alias{heterozygosity}
\alias{heterozygosity,SeqVarGDSClass-method}
\alias{homozygosity}
\alias{homozygosity,SeqVarGDSClass-method}
\alias{hethom}
\alias{hethom,SeqVarGDSClass-method}

\title{Heterozygosity and Homozygosity}

\description{
  Calculate heterozygosity and homozygosity by variant or by sample
}

\usage{
\S4method{heterozygosity}{SeqVarGDSClass}(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE,
    parallel=FALSE)
\S4method{homozygosity}{SeqVarGDSClass}(gdsobj, allele=c("any", "ref", "alt"), margin=c("by.variant", "by.sample"), use.names=FALSE,
    parallel=FALSE)
\S4method{hethom}{SeqVarGDSClass}(gdsobj, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{margin}{Possible values are "by.variant" or "by.sample,"
    indicating whether the calculation should be done over all samples
    for each variant, or over all variants for each sample.}
  \item{use.names}{A logical indicating whether to assign variant or
    samples IDs as names of the output vector.}
  \item{allele}{Possible values are "any", "ref," or "alt," indicating
    which alleles to consider when calculating homozygosity.}
  \item{parallel}{Logical, numeric, or other value to control parallel
    processing; see \code{\link{seqParallel}} for details. Only applies
  if \code{margin="as.variant"}.}
}

\details{
  \code{heterozyogosity} calulates the fraction of heterozygous genotypes
  in a GDS object, either by variant or by sample.

  \code{homozygosity} calculates the rate of homozygous genotypes in a
  GDS object, either by sample or by variant.  If \code{allele="any"},
  all homozygous genotypes are considered (reference or any alternate
  allele).  If \code{allele="ref"}, only reference homozygotes are
  considered.  If \code{allele="alt"}, any alternate allele homozygote
  is considered.  For example, "ref" will count "0/0" genotypes only,
  "alt" will count "1/1", "2/2", etc. (but not "0/0"), and "any" will
  count all of the above.

  \code{hethom} calculates the ratio of heterozygous genotypes to alternate homozygous genotypes by sample.
}

\value{
  A numeric vector of heterozyogity or homozygosity rates.  If
  \code{margin="by.variant"}, the vector will have length equal to the
  number of variants in the GDS object.  If \code{margin="by.sample"},
  the vector will have length equal to the number of samples.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{alleleFrequency}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
head(heterozygosity(gds, margin="by.variant"))
head(homozygosity(gds, allele="any", margin="by.variant"))
head(homozygosity(gds, allele="ref", margin="by.variant"))
head(homozygosity(gds, allele="alt", margin="by.variant"))

## Het/Hom Non-Ref by sample
head(hethom(gds))

seqClose(gds)
}

\keyword{manip}
