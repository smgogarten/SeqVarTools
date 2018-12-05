\name{imputedDosage}
\alias{imputedDosage}
\alias{imputedDosage,SeqVarGDSClass-method}

\title{Get imputed dosage}

\description{
  Get matrix of imputed dosage values from a GDS object
}

\usage{
\S4method{imputedDosage}{SeqVarGDSClass}(gdsobj, dosage.field="DS", use.names=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{dosage.field}{The name of the dosage field in the GDS object (will be prepended with \code{"annotation/format"}).}
  \item{use.names}{A logical indicating whether to assign sample and
    variant IDs as dimnames of the resulting matrix.}
}

\details{
  Reads dosage from the dosage-specific field in the GDS object, rather than counting alleles from called genotypes.

  Only one dosage value per variant is allowed; the method will return an error if multiple dosages are present for a single variant.
}

\value{
  A numeric matrix of dosage values with dimensions [sample,variant].
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{refDosage}},
  \code{\link{altDosage}}
}

\examples{
# convert VCF to GDS keeping dosage field
vcffile <- system.file("extdata", "gl_chr1.vcf", package="SeqVarTools")
gdsfile <- tempfile()
seqVCF2GDS(vcffile, gdsfile, fmt.import="DS", verbose=FALSE)

gds <- seqOpen(gdsfile)
dos <- imputedDosage(gds)
head(dos)
seqClose(gds)
unlink(gdsfile)
}

\keyword{manip}

