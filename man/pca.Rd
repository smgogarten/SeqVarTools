\name{pca}
\alias{pca}
\alias{pca,SeqVarGDSClass-method}

\title{Principal Component Analysis}

\description{Calculates the eigenvalues and eignevectors of a
  SeqVarGDSClass object with Principal Component Analysis}

\usage{
\S4method{pca}{SeqVarGDSClass}(gdsobj, eigen.cnt=32)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{eigen.cnt}{An integer indicating how many eigenvalues and
    eignvectors to return.}
}

\details{
  Calculates the genetic covariance matrix and finds the eigen decomposition.
}

\value{
  A list with two elements:
  \item{eigenval}{A vector of length \code{eigen.cnt} with eigenvalues}
  \item{eigenvect}{A matrix of dimension ("selected samples", \code{eigen.cnt}).}
}

\author{Xiuwen Zheng, Stephanie Gogarten}

\references{
  Patterson N, Price AL, Reich D (2006) Population structure and
  eigenanalysis. PLoS Genetics 2:e190.
}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
pca <- pca(gds)
pca$eigenval
head(pca$eigenvect)
seqClose(gds)
}

\keyword{manip}
