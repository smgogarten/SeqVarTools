\name{pedigree}
\alias{pedigree}
\docType{data}
\title{Pedigree for example data}
\description{
  Pedigree for example data files in SeqArray.
}
\usage{
pedigree
}
\format{
  A data.frame with the following columns.
  \describe{
    \item{\code{family}}{Family ID}
    \item{\code{individ}}{Individual ID}
    \item{\code{father}}{Father ID}
    \item{\code{mother}}{Mother ID}
    \item{\code{sex}}{Sex}
    \item{\code{sample.id}}{sample.id in VCF/GDS files}
  }
}
\details{
  There is one trio in the pedigree.
}
\source{
HapMap
}
\examples{
data(pedigree)
head(pedigree)
gds <- seqOpen(seqExampleFileName("gds"))
setdiff(seqGetData(gds, "sample.id"), pedigree$sample.id)
seqClose(gds)
}
\keyword{datasets}
