\name{meanBySample}
\alias{meanBySample}
\alias{meanBySample,SeqVarGDSClass-method}

\title{Mean value by sample}

\description{Calculate the mean value of a variable by sample over all variants}

\usage{
\S4method{meanBySample}{SeqVarGDSClass}(gdsobj, var.name, use.names=FALSE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{var.name}{Character string with name of the variable, most
    likely "annotation/format/VARIABLE_NAME".}
  \item{use.names}{A logical indicating whether to assign sample IDs
    as names of the output vector.}
}

\details{
  Mean values by variant can be calculated using \code{seqApply(gdsobj, var.name,
    mean, na.rm=TRUE)}.  Currently \code{seqApply} can only be used with
  the option \code{margin="by.variant"}.
  This method provides a way to calculate mean values by sample.
}

\value{
  A numeric vector of mean values.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{seqApply}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
head(meanBySample(gds, "annotation/format/DP", use.names=TRUE))
seqClose(gds)
}

\keyword{manip}
