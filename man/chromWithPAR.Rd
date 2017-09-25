\name{chromWithPAR}
\alias{chromWithPAR}
\alias{chromWithPAR,SeqVarGDSClass-method}

\title{Identify pseudoautosomal region}

\description{Flag single nucleotide variants}

\usage{
\S4method{chromWithPAR}{SeqVarGDSClass}(gdsobj, genome.build=c("hg19", "hg38"))
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{genome.build}{A character sting indicating genome build.}
}

\details{
  The pseudoautosomal region (PAR) should be treated like the autosomes for purposes of calculating allele frequency. This method returns a vector where sex chromosome variants are labeled wither "X", "Y", or "PAR".
}

\value{
  A character vector of chromosome, with values "PAR" for the pseudoautosomal region.
}

\references{
\url{https://www.ncbi.nlm.nih.gov/grc/human}
}

\author{Stephanie Gogarten}
