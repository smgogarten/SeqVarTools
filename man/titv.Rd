\name{titv}
\alias{titv}
\alias{titv,SeqVarGDSClass-method}

\title{Transition/Transversion Ratio}

\description{
  Calculate transition/transversion ratio overall or by sample
}

\usage{
\S4method{titv}{SeqVarGDSClass}(gdsobj, by.sample=FALSE, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{by.sample}{A logical indicating whether TiTv should be
    calculated by sample or overall for the entire GDS object.}
  \item{use.names}{A logical indicating whether to assign
    sample IDs as names of the output vector (if \code{by.sample=TRUE}).}
}

\details{
  If \code{by.sample=FALSE} (the default), \code{titv} calulates the
  transition/transversion ratio (TiTv) over all samples.

  If \code{by.sample=TRUE}, \code{titv} calculates TiTv over all
  variant genotypes (heterozygous or homozygous non-reference) for each
  sample.
}

\value{
  A single value for TiTv if \code{by.sample=FALSE}.  If \code{by.sample=TRUE},
  a numeric vector containing TiTv for each sample.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{isVariant}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
titv(gds)
titv(gds, by.sample=TRUE)

## apply to a subset of variants
library(GenomicRanges)
chrom <- seqGetData(gds, "chromosome")
pos22 <- seqGetData(gds, "position")[chrom == 22]
ranges <- GRanges(seqnames="22", IRanges(min(pos22), max(pos22)))
applyMethod(gds, titv, ranges)

seqClose(gds)
}

\keyword{manip}
