\name{inbreedCoeff}
\alias{inbreedCoeff}
\alias{inbreedCoeff,SeqVarGDSClass-method}

\title{Inbreeding coefficient}

\description{Calculates the inbreeding coefficient by variant or by sample}

\usage{
\S4method{inbreedCoeff}{SeqVarGDSClass}(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE,
    parallel=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{margin}{Possible values are "by.variant" or "by.sample,"
    indicating whether the calculation should be done over all samples
    for each variant, or over all variants for each sample.}
  \item{use.names}{A logical indicating whether to assign variant or
    sample IDs as names of the output vector.}
  \item{parallel}{Logical, numeric, or other value to control parallel
    processing; see \code{\link{seqParallel}} for details. Only applies
  if \code{margin="as.variant"}.}
}

\details{
  For inbreeding coefficients by variant, calculates 1 - observed
  heterozygosity / expected heterozygosity.
  
  For individual inbreeding coefficients (\code{margin="by.sample"}),
  calculates Visscher's estimator described in Yang et al. (2010).
}

\value{
  Values for the inbreeding coefficient.
}

\author{Xiuwen Zheng, Stephanie Gogarten}

\references{
  Yang J, Benyamin B, McEvoy BP, Gordon S, Henders AK, Nyholt DR, Madden
  PA, Heath AC, Martin NG, Montgomery GW, Goddard ME, Visscher
  PM. 2010. Common SNPs explain a large proportion of the heritability
  for human height. Nat Genet. 42(7):565-9. Epub 2010 Jun 20.
}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
f <- inbreedCoeff(gds, margin="by.variant")
range(f, na.rm=TRUE)

ic <- inbreedCoeff(gds, margin="by.sample")
range(ic)
seqClose(gds)
}

\keyword{manip}
