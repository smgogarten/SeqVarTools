\name{getGenotype}
\alias{getGenotype}
\alias{getGenotypeAlleles}
\alias{refDosage}
\alias{altDosage}
\alias{alleleDosage}
\alias{getGenotype,SeqVarGDSClass-method}
\alias{getGenotypeAlleles,SeqVarGDSClass-method}
\alias{refDosage,SeqVarGDSClass-method}
\alias{altDosage,SeqVarGDSClass-method}
\alias{alleleDosage,SeqVarGDSClass,numeric-method}
\alias{alleleDosage,SeqVarGDSClass,list-method}

\title{Get genotype data}

\description{
  Get matrix of genotype values from a GDS object as VCF-style character
  strings
}

\usage{
\S4method{getGenotype}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
\S4method{getGenotypeAlleles}{SeqVarGDSClass}(gdsobj, use.names=TRUE, sort=FALSE)
\S4method{refDosage}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
\S4method{altDosage}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
\S4method{alleleDosage}{SeqVarGDSClass,numeric}(gdsobj, n=0, use.names=TRUE)
\S4method{alleleDosage}{SeqVarGDSClass,list}(gdsobj, n, use.names=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{use.names}{A logical indicating whether to assign sample and
    variant IDs as dimnames of the resulting matrix.}
  \item{sort}{Logical for whether to sort alleles lexographically ("G/T"
    instead of "T/G").}
  \item{n}{An integer, vector, or list indicating which allele(s) to return dosage
  for.  \code{n=0} is the reference allele, \code{n=1} is the first
  alternate allele, and so on.}
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

  If the argument \code{n} to\code{alleleDosage} is a single integer, the same allele is counted for all variants. If \code{n} is a vector with length=number of variants in the current filter, a different allele is counted for each variant. If \code{n} is a list, more than one allele can be counted for each variant. For example, if \code{n[[1]]=c(1,3)}, genotypes "0/1" and "0/3" will each have a dosage of 1 and genotype "1/3" will have a dosage of 2.
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

  \code{altDosage} returns an integer matrix with the dosage of any
  alternate allele: 2 for two alternate alleles ("1/1", "1/2", etc.), 1
  for one alternate allele, and 0 for no alternate allele (homozygous reference).

  \code{alleleDosage} with an integer argument returns an integer matrix with the dosage of the
  specified allele only: 2 for two copies of the allele ("0/0" if \code{n=0}, "1/1" if \code{n=1}, etc.), 1
  for one copy of the specified allele, and 0 for no copies of the allele.

  \code{alleleDosage} with a list argument returns a list of allele x sample matrices with the dosage of each specified allele for each variant.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link[SeqArray]{seqGetData}},
  \code{\link[SeqArray]{seqSetFilter}},
  \code{\link{alleleFrequency}}
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
altDosage(gds)
alleleDosage(gds, n=0)
alleleDosage(gds, n=1)
alleleDosage(gds, n=c(0,1,0,1,0))
alleleDosage(gds, n=list(0,c(0,1),0,c(0,1),1))
seqClose(gds)
}

\keyword{manip}

