\name{hwe}
\alias{hwe}
\alias{hwe,SeqVarGDSClass-method}

\title{Exact test for Hardy-Weinberg equilibrium}

\description{Performs an exact test for Hardy-Weinberg equilibrium on
  Single-Nucleotide Variants}

\usage{
\S4method{hwe}{SeqVarGDSClass}(gdsobj, permute=FALSE, parallel=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{permute}{A logical indicating whether to permute the genotypes to get a set of p-values under the null hypothesis.}
  \item{parallel}{Logical, numeric, or other value to control parallel
    processing; see \code{\link{seqParallel}} for details.}
}

\details{
  HWE calculations are performed with the \code{\link{HWExact}}
  function in the \pkg{\link{GWASExactHW}} package.

  \code{permute=TRUE} will permute the genotypes prior to running the test. This can be useful for obtaining a set of expected values under the null hypothesis to compare to the observed values.

  P values are set to \code{NA} for all multiallelic
  and monomorphic variants.
}

\value{
  A data.frame with the following columns:
\item{variant.id}{The unique identifier for the variant}
\item{nAA}{The number of reference homozygotes}
\item{nAa}{The number of heterozygotes}
\item{naa}{The number of alternate homozygotes}
\item{afreq}{The reference allele frequency}
\item{p}{p values for the exact test}
\item{f}{The inbreeding coefficient, 1 - observed
  heterozygosity / expected heterozygosity}
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
hw <- applyMethod(gds, hwe, variant=var.auto)
head(hw)
sum(is.na(hw$p))
range(hw$p, na.rm=TRUE)
seqClose(gds)
}

\keyword{manip}
