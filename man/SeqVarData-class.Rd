\name{SeqVarData}
\docType{class}

\alias{SeqVarData-class}
\alias{SeqVarData}
\alias{show,SeqVarData-method}
\alias{sampleData}
\alias{sampleData,SeqVarData-method}
\alias{sampleData<-}
\alias{sampleData<-,SeqVarData,AnnotatedDataFrame-method}
\alias{variantData}
\alias{variantData,SeqVarData-method}
\alias{variantData<-}
\alias{variantData<-,SeqVarData,AnnotatedDataFrame-method}
\alias{granges}
\alias{granges,SeqVarData-method}
\alias{validateSex}
\alias{validateSex,SeqVarData-method}

\title{SeqVarData}

\description{
	Extends \code{\link{SeqVarGDSClass}} to include annotation for samples and variants.
}

\details{
  A \code{SeqVarData} object adds an \code{\link{AnnotatedDataFrame}} for both samples and variants to a \code{\link{SeqVarGDSClass}} object.

  Note that a \code{SeqVarData} object must be created using an \emph{unfiltered}  \code{\link{SeqVarGDSClass}} object. The \code{sample.id} column in the \code{sampleData} \code{\link{AnnotatedDataFrame}} must exactly match the \code{sample.id} node in the GDS file (and similarly for \code{variant.id} in \code{variantData}). This enables all subsequent filters set on the \code{SeqVarData} object to apply to the GDS and the annotation simultaneously.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{SeqVarData(gds, sampleData, variantData)}: Returns a \code{SeqVarData} object.

      \code{gds} can be either the filename of a sequencing GDS file or an existing \code{\link{SeqVarGDSClass}} object.

      \code{sampleData} must be an \code{\link{AnnotatedDataFrame}} with a column \code{sample.id} matching \code{sample.id} in the GDS file. If this argument is missing, a data frame with 0 columns will be created.

      \code{variantData} must be an \code{\link{AnnotatedDataFrame}} with a column \code{variant.id} matching \code{variant.id} in the GDS file. If this argument is missing, a data frame with 0 columns will be created.
    }
  }
}


\section{Accessors}{
	\describe{	
		\item{}{
			\code{sampleData(x)}, \code{sampleData(x)<- value}:
				Get or set the \code{\link{AnnotatedDataFrame}} with sample data. If a sample filter has been applied with \code{\link{seqSetFilter}}, only selected samples will be returned. \code{value} must include all samples.
		}	
		\item{}{
			\code{variantData(x)}, \code{variantData(x)<- value}:
				Get or set the \code{\link{AnnotatedDataFrame}} with variant data. If a variant filter has been applied with \code{\link{seqSetFilter}}, only selected variants will be returned. \code{value} must include all variants.
		}	
		\item{}{
			\code{granges(x)}:
				Return a \code{\link{GRanges}} object with the columns of \code{variantData} as metadata columns.
		}	
		\item{}{
			\code{validateSex(x)}:
				Return the contents of a column named "sex" in \code{sampleData(x)}, provided the contents are valid (values either "M"/"F" or 1/2, or NA). If the column is missing or invalid, return \code{NULL}.
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

head(validateSex(seqData))

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


