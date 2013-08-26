\name{getGenotype}
\alias{getGenotype}
\alias{getGenotypeAlleles}
\alias{refDosage}
\alias{getGenotype,SeqVarGDSClass-method}
\alias{getGenotypeAlleles,SeqVarGDSClass-method}
\alias{refDosage,SeqVarGDSClass-method}

\title{Get genotype data}

\description{
  Get matrix of genotype values from a GDS object as VCF-style character
  strings
}

\usage{
\S4method{getGenotype}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
\S4method{getGenotypeAlleles}{SeqVarGDSClass}(gdsobj, use.names=TRUE, sort=FALSE)
\S4method{refDosage}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{use.names}{A logical indicating whether to assign sample and
    variant IDs as dimnames of the resulting matrix.}
  \item{sort}{Logical for whether to sort alleles lexographically ("G/T"
    instead of "T/G").}
}

\details{
  In \code{getGenotype}, genotypes are coded as in the VCF file, where "0/0" is homozygous
  reference, "0/1" is heterozygous for the first alternate allele, "0/2"
  is heterozygous for the second alternate allele, etc.

  Separators are
  "/" for unphased and "|" for phased.  If \code{sort=TRUE}, all
  returned genotypes will be unphased.
  Missing genotypes are coded as \code{NA}.
  
  Only diploid genotypes (the first two alleles at a given site) are returned.
}

\value{
  \code{getGenotype} and \code{getGenotypeAlleles} return a character matrix with dimensions [sample,variant] containing diploid
  genotypes.

  \code{getGenotype} returns alleles as "0", "1", "2", etc. indicating
  refernence and alternate alleles.

  \code{getGenotypeAlleles} returns alleles as "A", "C", "G", "T".
  \code{sort=TRUE} sorts lexographically, which may be useful for
  comparing genotypes with data generated using a different reference
  sequence.

  \code{refDosage} returns an integer matrix with the dosage of the
  reference allele: 2 for two copies of the reference allele ("0/0"), 1
  for one copy of the reference allele, and 0 for two alternate alleles.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link[SeqArray]{seqGetData}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
variant.id <- seqGetData(gds, "variant.id")
sample.id <- seqGetData(gds, "sample.id")
seqSetFilter(gds, variant.id=variant.id[1:5],
             sample.id=sample.id[1:10])
getGenotype(gds)
getGenotypeAlleles(gds)
refDosage(gds)
seqClose(gds)
}

\keyword{manip}

