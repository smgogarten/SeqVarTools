\name{countSingletons}
\alias{countSingletons}
\alias{countSingletons,SeqVarGDSClass-method}

\title{Count singletons}

\description{Count singleton variants for each sample}

\usage{
\S4method{countSingletons}{SeqVarGDSClass}(gdsobj, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
}

\details{
  A singleton variant is a variant in which only one sample has a non-reference allele.  For each sample, \code{countSingletons} finds the number of variants for which that sample has the only non-reference allele.
}

\value{
  A vector of the number of singleton variants per sample.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{alleleFrequency}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
head(countSingletons(gds))
seqClose(gds)
}

\keyword{manip}
