\name{alleleFrequency}
\alias{alleleFrequency}
\alias{alleleFrequency,SeqVarGDSClass-method}
\alias{alleleFrequency,SeqVarData-method}
\alias{alleleCount}
\alias{alleleCount,SeqVarGDSClass-method}
\alias{alleleCount,SeqVarData-method}
\alias{minorAlleleCount}
\alias{minorAlleleCount,SeqVarGDSClass-method}
\alias{minorAlleleCount,SeqVarData-method}

\title{Allele frequency}

\description{Calculate allele frequency for each variant}

\usage{
\S4method{alleleFrequency}{SeqVarGDSClass}(gdsobj, n=0, use.names=FALSE, parallel=FALSE)
\S4method{alleleFrequency}{SeqVarData}(gdsobj, n=0, use.names=FALSE, sex.adjust=TRUE, male.diploid=TRUE,
    genome.build=c("hg19", "hg38"), parallel=FALSE)
\S4method{alleleCount}{SeqVarGDSClass}(gdsobj, n=0, use.names=FALSE, parallel=FALSE)
\S4method{alleleCount}{SeqVarData}(gdsobj, n=0, use.names=FALSE, sex.adjust=TRUE, male.diploid=TRUE,
    genome.build=c("hg19", "hg38"), parallel=FALSE)
\S4method{minorAlleleCount}{SeqVarData}(gdsobj, use.names=FALSE, sex.adjust=TRUE, male.diploid=TRUE,
    genome.build=c("hg19", "hg38"), parallel=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{n}{An integer indicating which allele to calculate the frequency
  of.  \code{n=0} is the reference allele, \code{n=1} is the first
  alternate allele, and so on.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
  \item{sex.adjust}{Logical for whether to adjust frequency calculations based on sex. If \code{TRUE}, X chromosome frequency (excluding the PAR) will be calculated assuming the dosage of the specifed allele for males is half that for females. Y chromosome frequency will be calculated using males only.}
  \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid.}
  \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes.}
  \item{parallel}{Logical, numeric, or other value to control parallel
    processing; see \code{\link{seqParallel}} for details.}
}

\details{
  Frequency or count can be calculated over any allele, specified by the argument
  \code{n}.  Default is the reference allele frequency (\code{n=0}).

  The \code{\link{SeqVarData}} method will calculate frequency and count correctly for X and Y chromosomes, provided a column "sex" is included in the sampleData slot with values "M"/"F" or 1/2. Arguments given to this method are passed to the parent method for \code{\link{SeqVarGDSClass}}. If the ploidy of the "genotype" node in the GDS file is 2, the default assumption is that genotypes for males on sex chromosomes are coded as diploid, "0/0" or "1/1". If this is not the case, use \code{male.diploid=FALSE}.

  For multiallelic variants, the minor allele count will be the smaller of the reference allele count or the sum of all alternate allele counts. 
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
