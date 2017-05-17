\name{isSNV}
\alias{isSNV}
\alias{isSNV,SeqVarGDSClass-method}

\title{Flag single nucleotide variants}

\description{Flag single nucleotide variants}

\usage{
\S4method{isSNV}{SeqVarGDSClass}(gdsobj, biallelic=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{biallelic}{A logical indicating whether only biallelic SNVs are considered.}
}

\details{
  If \code{biallelic=TRUE}, a variant is considered a single nucleotide variant (SNV) if there is
  one reference allele and one alternate allele, each one base in
  length.  If \code{biallelic=FALSE}, there may be multiple alternate
  alleles, each one base in length.

  Setting \code{biallelic=TRUE} is considerably faster for large data sets.
}

\value{
  A logical vector indicating which variants are SNVs.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{allele-methods}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
table(isSNV(gds))
seqClose(gds)
}
