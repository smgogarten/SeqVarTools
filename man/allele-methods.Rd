\name{allele-methods}
\alias{allele-methods}
\alias{refChar}
\alias{refChar,SeqVarGDSClass-method}
\alias{altChar}
\alias{altChar,SeqVarGDSClass-method}
\alias{nAlleles}
\alias{nAlleles,SeqVarGDSClass-method}

\title{Extract allele information from a GDS object}

\description{
  Extract reference and alternate alleles and allele counts
  from a GDS object.
}

\usage{
\S4method{refChar}{SeqVarGDSClass}(gdsobj)
\S4method{altChar}{SeqVarGDSClass}(gdsobj, n=0)
\S4method{nAlleles}{SeqVarGDSClass}(gdsobj)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{n}{An integer indicating which alternate allele to return.  \code{n=0} returns a
    comma-separated string of all alternate alleles.}
}

\details{
  These methods parse the "allele" field of a GDS object.
}

\value{
  \code{refChar} returns a character vector of reference alleles.

  \code{altChar} returns a character vector of alternate alleles.  If
  \code{n=0}, multiple alternate alleles are represented as a
  comma-separated string.  If \code{n>0}, only the \code{n}th alternate
  allele is returned.

  \code{nAlleles} returns an integer vector of the number of alleles
  (reference and alternate) for each variant.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
table(refChar(gds))
table(altChar(gds))
table(altChar(gds, n=1))
table(altChar(gds, n=2), useNA="ifany")
table(nAlleles(gds))
seqClose(gds)
}

\keyword{manip}
