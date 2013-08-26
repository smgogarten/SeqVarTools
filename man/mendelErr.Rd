\name{mendelErr}
\alias{mendelErr}
\alias{mendelErr,SeqVarGDSClass-method}

\title{Mendelian errors}

\description{Detect Mendelian errors}

\usage{
\S4method{mendelErr}{SeqVarGDSClass}(gdsobj, pedigree, use.names=FALSE,
autosomes=1:22, xchrom="X", ychrom="Y", verbose=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{pedigree}{A data.frame with columns (family, individ, father,
  mother, sex, sample.id).  "sex" column should have values "M"/"F".
  "sample.id" values should correspond to "sample.id" in \code{gdsobj}.}
  \item{use.names}{A logical indicating whether to assign variant IDs
    as names of the output vector.}
  \item{autosomes}{A vector with chromosome values in \code{gdsobj} corresponding to autosomes.}
  \item{xchrom}{The chromosome value in \code{gdsobj} corresponding to
    the X chromosome.}
  \item{ychrom}{The chromosome value in \code{gdsobj} corresponding to
    the Y chromosome.}
  \item{verbose}{A logical indicating whether to print the number of
  samples selected for each trio.}
}

\details{
  Mendelian errors are detected for each trio in \code{pedigree}.  Duos
  (mother or father missing) are included.  The pedigree must have only
  one sample per individual.
}

\value{
  A list with the following elements:
  \item{by.variant}{An integer vector with the number of mendelian errors detected for
    each variant.  If \code{use.names=TRUE}, the vector will be named
  with variant IDs.}
  \item{by.trio}{An integer vector with the number of mendelian errors
  detected for each trio.  The vector will be named with the sample ID
  of the child in each trio.}
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
data(pedigree)
err <- mendelErr(gds, pedigree)
table(err$by.variant)
err$by.trio
seqClose(gds)
}

\keyword{manip}
