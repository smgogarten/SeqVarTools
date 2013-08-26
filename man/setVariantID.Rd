\name{setVariantID}
\alias{setVariantID}

\title{Change the variant ID of a GDS file}

\description{
  Replace the variable "variant.id" in a GDS file with a user-supplied
  unique vector of the same length.
}

\usage{
setVariantID(gdsfile, variant.id)
}

\arguments{
  \item{gdsfile}{A character string with the file path of a GDS file.}
  \item{variant.id}{A vector with the new variant IDs.}
}

\details{
  A VCF file created by \code{\link{seqVCF2GDS}} creates a variable
  "variant.id" containing sequential integers to identify each variant.
  \code{setVariantID} allows the user to replace these values with
  something more meaningful.  The replacement values in
  \code{variant.id} must be unique and have the same length as the
  original "variant.id" vector.

  Using character values for \code{variant.id} may affect performance
  for large datasets.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{seqVCF2GDS}}
}

\examples{
oldfile <- system.file("extdata", "gl_chr1.gds", package="SeqVarTools")
newfile <- tempfile()
file.copy(oldfile, newfile)

gds <- seqOpen(newfile)
rsID <- seqGetData(gds, "annotation/id")
seqClose(gds)

setVariantID(newfile, rsID)
gds <- seqOpen(newfile)
seqGetData(gds, "variant.id")
head(getGenotype(gds))
seqClose(gds)
 
unlink(newfile)
}

\keyword{manip}
