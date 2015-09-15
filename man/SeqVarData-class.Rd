\name{SeqVarData}
\docType{class}

\alias{SeqVarData-class}
\alias{SeqVarData}
\alias{sampleData}
\alias{sampleData,SeqVarData-method}
\alias{sampleData<-}
\alias{sampleData<-,SeqVarData,AnnotatedDataFrame-method}

\title{SeqVarData}

\description{
	Extends \code{\link{SeqVarGDSClass}} to include annotation for the samples.
}

\details{
	A \code{SeqVarData} object adds an \code{\link{AnnotatedDataFrame}} to a \code{\link{SeqVarGDSClass}} object.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{SeqVarData(gds, sampleData)}: Returns a \code{SeqVarData} object.

      \code{gds} can be either the filename of a sequencing GDS file or an existing \code{\link{SeqVarGDSClass}} object.

      \code{sampleData} must be an \code{\link{AnnotatedDataFrame}} with a column \code{sample.id} matching \code{sample.id} in the GDS file.
    }
  }
}


\section{Accessors}{
	\describe{	
		\item{}{
			\code{sampleData(x)}, \code{sampleData(x)<- value}:
				Get or set the \code{\link{AnnotatedDataFrame}} with sample data. If a sample filter has been applied with \code{\link{seqSetFilter}}, only selected samples will be returned. \code{value} must include all samples.
		}
	}

	See \code{\link{SeqVarGDSClass}} for additional methods.
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{seqVCF2GDS}},
  \code{\link{seqOpen}},
  \code{\link{seqGetData}},
  \code{\link{seqSetFilter}},
  \code{\link{seqApply}},
  \code{\link{seqClose}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))

## create sample annotation
library(Biobase)
sample.id <- seqGetData(gds, "sample.id")
sex <- sample(c("M","F"), length(sample.id), replace=TRUE)
phenotype <- rnorm(length(sample.id), mean=10)
samp <- data.frame(sample.id, sex, phenotype, stringsAsFactors=FALSE)
meta <- data.frame(labelDescription=c("unique sample identifier",
     "sex (M=male, f=female)", "example phenotype"), 
      row.names=names(samp), stringsAsFactors=FALSE)
sample.data <- AnnotatedDataFrame(samp, meta)

seqData <- SeqVarData(gds, sample.data)

## add another annotation column
sample.data$site <- sample(letters, length(sample.id), replace=TRUE)
varMetadata(sample.data)["site", "labelDescription"] <- "study site"
sampleData(seqData) <- sample.data

## set a filter
seqSetFilter(seqData, sample.id=sample.id[1:10])
nrow(sampleData(seqData))

seqClose(seqData)
}

\keyword{methods}
\keyword{classes}


