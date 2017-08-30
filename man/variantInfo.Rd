\name{variantInfo}
\alias{variantInfo}
\alias{variantInfo,SeqVarGDSClass-method}
\alias{expandedVariantIndex}
\alias{expandedVariantIndex,SeqVarGDSClass-method}

\title{Variant info}

\description{Return basic variant info as a data.frame.}

\usage{
\S4method{variantInfo}{SeqVarGDSClass}(gdsobj, alleles=TRUE, expanded=FALSE)
\S4method{expandedVariantIndex}{SeqVarGDSClass}(gdsobj)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{alleles}{A logical value for whether to include ref and alt alleles}
  \item{expanded}{A logical value for whether to expand multi-allelic variants with one row per alternate allele.}
}

\details{
Variants can be represented in collapsed form, with one row per variant, or in expanded form, with one row per alternate allele for multiallelic variants.
}

\value{
\code{variantInfo} returns a data.frame with variant.id, chromosome, and position for each variant. If \code{alleles=TRUE}, the data.frame includes ref and alt. If \code{expanded=TRUE}, the data.frame includes allele.index, which is 1 for the first alternate allele, 2 for the second alternate, etc.

\code{expandedVariantIndex} returns an index to transform a vector or matrix from collapsed to expanded form.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
seqSetFilter(gds, variant.sel=1323:1327)
variantInfo(gds, alleles=TRUE)
variantInfo(gds, alleles=TRUE, expanded=TRUE)
expandedVariantIndex(gds)
seqClose(gds)
}
