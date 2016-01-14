\name{regression}
\alias{regression}
\alias{regression,SeqVarData-method}

\title{Linear or logistic regression}

\description{Run linear or logistic regression on variants}

\usage{
\S4method{regression}{SeqVarData}(gdsobj, outcome, covar=NULL,
                                  model.type=c("linear", "logistic", "firth"))
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarData}} object}
  \item{outcome}{A character string with the name of the column in \code{sampleData(gdsobj)} containing the outcome variable}
  \item{covar}{A character vector with the name of the column(s) in \code{sampleData(gdsobj)} containing the covariates}
  \item{model.type}{the type of model to be run.  "linear" uses \code{\link{lm}}, "logistic" uses \code{\link{glm}} with \code{family=binomial()}, and "firth" uses \code{\link{logistf}}.}
}

%\details{
%}

\value{
a data.frame with the following columns (if applicable):
  \item{variant.id}{variant identifier}
  \item{n}{number of samples with non-missing data}
  \item{n0}{number of controls (outcome=0) with non-missing data}
  \item{n1}{number of cases (outcome=1) with non-missing data}
  \item{freq}{reference allele frequency}
  \item{freq0}{reference allele frequency in controls}
  \item{freq1}{reference allele frequency in cases}
  \item{Est}{beta estimate for genotype}
  \item{SE}{standard error of beta estimate for the genotype}
  \item{Wald.Stat}{chi-squared test statistic for association}
  \item{Wald.pval}{p-value for association}
  \item{PPL.Stat}{firth only: profile penalized likelihood test statistic for association}
  \item{PPL.pval}{firth only: p-value for association}
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarData}},
  \code{\link{seqSetFilter}},
  \code{\link{lm}},
  \code{\link{glm}},
  \code{\link{logistf}}	
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))

## create some phenotype data
library(Biobase)
sample.id <- seqGetData(gds, "sample.id")
n <- length(sample.id)
df <- data.frame(sample.id,
   sex=sample(c("M", "F"), n, replace=TRUE),
   age=sample(18:70, n, replace=TRUE),
   phen=rnorm(n),
   stringsAsFactors=FALSE)
meta <- data.frame(labelDescription=c("sample identifier",
   "sex", "age", "phenotype"), row.names=names(df))
sample.data <- AnnotatedDataFrame(df, meta)
seqData <- SeqVarData(gds, sample.data)

## select samples and variants
seqSetFilter(gds, sample.id=sample.id[1:50], variant.id=1:10)

res <- regression(seqData, outcome="phen", covar=c("sex", "age"))
res
seqClose(gds)
}

\keyword{manip}
