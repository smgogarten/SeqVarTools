\name{isVariant}
\alias{isVariant}
\alias{isVariant,SeqVarGDSClass-method}

\title{Locate variant samples across sites}

\description{
  Locate which samples are variant for each site in a GDS object
}

\usage{
\S4method{isVariant}{SeqVarGDSClass}(gdsobj, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{use.names}{A logical indicating whether to assign sample and
    variant IDs as dimnames of the resulting matrix.}
}

\details{
  Each sample/site cell of the resulting matrix is \code{TRUE} if the
  genotype at that location for that sample contains an alternate
  allele.  A genotype of "0/0" is not variant, while genotypes
  "0/1", "1/0", "0/2", etc. are variant.
}

\value{
  A logical matrix with dimensions [sample,site] which is \code{TRUE}
  for cells where the genotype contains an alternate allele.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{getGenotype}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
variant.id <- seqGetData(gds, "variant.id")
sample.id <- seqGetData(gds, "sample.id")
applyMethod(gds, isVariant, variant.id[1:5], sample.id[1:10])
applyMethod(gds, isVariant, variant.id[1:5], sample.id[1:10], use.names=TRUE)
seqClose(gds)
}

\keyword{manip}
