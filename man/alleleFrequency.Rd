\name{alleleFrequency}
\alias{alleleFrequency}
\alias{alleleFrequency,SeqVarGDSClass-method}

\title{Allele frequency}

\description{Calculate allele frequency for each variant}

\usage{
\S4method{alleleFrequency}{SeqVarGDSClass}(gdsobj, n=0, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{n}{An integer indicating which allele to calculate the frequency
  of.  \code{n=0} is the reference allele, \code{n=1} is the first
  alternate allele, and so on.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
}

\details{
  Frequency can be calculated over any allele, specified by the argument
  \code{n}.  Default is the reference allele frequency (\code{n=0}).
}

\value{
  A numeric vector of allele frequencies.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{heterozygosity}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
head(alleleFrequency(gds))
head(alleleFrequency(gds, n=1))
head(alleleFrequency(gds, n=2))
seqClose(gds)
}

\keyword{manip}
