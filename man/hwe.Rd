\name{hwe}
\alias{hwe}
\alias{hwe,SeqVarGDSClass-method}

\title{Exact test for Hardy-Weinberg equilibrium}

\description{Performs an exact test for Hardy-Weinberg equilibrium on
  Single-Nucleotide Variants}

\usage{
\S4method{hwe}{SeqVarGDSClass}(gdsobj, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
}

\details{
  HWE calculations are performed with the \code{\link{HWExact}}
  function in the \pkg{\link{GWASExactHW}} package.

  P values are set to \code{NA} for all non-single-nucleotide variants
  and monomorphic variants.
}

\value{
  A vector of p values for the exact test.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
## autosomal variants only
auto <- seqGetData(gds, "chromosome") \%in\% 1:22
var.auto <- seqGetData(gds, "variant.id")[auto]
pv <- applyMethod(gds, hwe, variant=var.auto)
head(pv)
sum(is.na(pv))
range(pv, na.rm=TRUE)
seqClose(gds)
}

\keyword{manip}
