\name{applyMethod}
\alias{applyMethod}
\alias{applyMethod,SeqVarGDSClass,function,character-method}
\alias{applyMethod,SeqVarGDSClass,function,numeric-method}
\alias{applyMethod,SeqVarGDSClass,function,GRanges-method}
\alias{applyMethod,SeqVarGDSClass,function,missing-method}

\title{Apply method to GDS object}

\description{Apply a method to a subset of variants and/or samples in a GDS object}

\usage{
\S4method{applyMethod}{SeqVarGDSClass,function,character}(gdsobj, FUN, variant, sample=NULL, ...)
\S4method{applyMethod}{SeqVarGDSClass,function,numeric}(gdsobj, FUN, variant, sample=NULL, ...)
\S4method{applyMethod}{SeqVarGDSClass,function,GRanges}(gdsobj, FUN, variant, sample=NULL, ...)
\S4method{applyMethod}{SeqVarGDSClass,function,missing}(gdsobj, FUN, variant, sample=NULL, ...)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{FUN}{A method or function to be applied to \code{gdsobj}.}
  \item{variant}{A vector of variant.id values or a \code{GRanges}
    object defining the variants to be included in the call to
    \code{FUN}.}
  \item{sample}{A vector of sample.id values defining the samples to be
    included in the call to \code{FUN}.}
  \item{\dots}{Additional arguments, passed to \code{FUN}.}
}

\details{
  \code{applyMethod} applies a method or function \code{FUN} to the
  subset of variants defined by \code{variant} and samples defined by
  \code{sample} in a GDS object.

  If a filter was previously set with \code{\link{seqSetFilter}}, it
  will be saved and reset after the call to \code{applyMethod}.
}

\value{
  The result of the call to \code{FUN}.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
variant.id <- seqGetData(gds, "variant.id")
sample.id <- seqGetData(gds, "sample.id")
applyMethod(gds, getGenotype, variant.id[1:5], sample.id[1:10])

library(GenomicRanges)
chrom <- seqGetData(gds, "chromosome")
pos22 <- seqGetData(gds, "position")[chrom == 22]
ranges <- GRanges(seqnames="22", IRanges(min(pos22), max(pos22)))
applyMethod(gds, heterozygosity, ranges, margin="by.sample")
applyMethod(gds, heterozygosity, ranges, margin="by.variant")

seqClose(gds)
}

\keyword{manip}
