\name{alleleFrequency}
\alias{alleleFrequency}
\alias{alleleFrequency,SeqVarGDSClass-method}
\alias{alleleFrequency,SeqVarData-method}

\title{Allele frequency}

\description{Calculate allele frequency for each variant}

\usage{
\S4method{alleleFrequency}{SeqVarGDSClass}(gdsobj, n=0, use.names=FALSE)
\S4method{alleleFrequency}{SeqVarData}(gdsobj, n=0, use.names=FALSE, sex.adjust=TRUE, genome.build=c("hg19", "hg38"))
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{n}{An integer indicating which allele to calculate the frequency
  of.  \code{n=0} is the reference allele, \code{n=1} is the first
  alternate allele, and so on.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
  \item{sex.adjust}{Logical for whether to adjust frequency calculations based on sex. If \code{TRUE}, X chromosome frequency (excluding the PAR) will be calculated assuming the dosage of the specifed allele for males is half that for females. Y chromosome frequency will be calculated using males only.}
  \item{genome.build}{A character sting indicating genome build.}
}

\details{
  Frequency can be calculated over any allele, specified by the argument
  \code{n}.  Default is the reference allele frequency (\code{n=0}).

  The \code{\link{SeqVarData}} method will calculate frequency correctly for X and Y chromosomes, provided a column "sex" is included in the sampleData slot with values "M"/"F" or 1/2. Arguments given to this method are passed to the parent method for \code{\link{SeqVarGDSClass}}.
}

\value{
  A numeric vector of allele frequencies.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{chromWithPAR}},
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
