\name{refFrac}
\alias{refFrac}
\alias{refFrac,SeqVarGDSClass-method}
\alias{refFracOverHets}
\alias{refFracOverHets,SeqVarGDSClass-method}
\alias{refFracPlot}
\alias{refFracPlot,SeqVarGDSClass-method}

\title{Reference allele fraction}

\description{
  Calculate fraction of reference allele reads
}

\usage{
\S4method{refFrac}{SeqVarGDSClass}(gdsobj, use.names=TRUE)
\S4method{refFracOverHets}{SeqVarGDSClass}(gdsobj, FUN=mean, use.names=TRUE)
\S4method{refFracPlot}{SeqVarGDSClass}(gdsobj, variant.id, highlight=NULL, ...)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{FUN}{The function to apply over heterozygote calls (mean or median).}
  \item{use.names}{A logical indicating whether to assign variant or
    samples IDs as names of the output vector.}
  \item{variant.id}{A vector of variant.ids to plot.}
  \item{highlight}{A list of sample.ids to highlight with sequential integers on each plot}
  \item{...}{Additional arguments passed to \code{\link{plot}}.}
}

\details{
  The variable "annotation/format/AD" (allelic depth) is required to
  compute the reference allele fraction.

  \code{refFracPlot} generates plots of total unfiltered depth (sum over "AD" for
  all alleles) versus reference allele fraction.  Points are color-coded
  by called genotype: teal = reference homozygote, orange = heterozygote
  including the reference allele, fuschia = heterozygote with two
  alternate alleles, purple = alternate homozygote, black = missing. 
  Darker colors indicate a higher density of points.
  Vertical black line
  is at 0.5, vertical teal line is the median reference allele fraction
  for ref/alt heterozygotes.  Values significantly different from 0.5
  (after applying a Bonferroni correction) are plotted with triangles.
}

\value{
  \code{refFrac} returns a sample by variant array of the reference allele
  fraction, defined as ref_depth / total_depth.

  \code{refFracOverHets} returns the mean (or other function,
  e.g. median) of reference allele
  depth (per variant) over all samples called as heterozygotes. 
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}}
}

\examples{
gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
gds <- seqOpen(gdsfile)
RF <- refFrac(gds)
dim(RF)
samples <- seqGetData(gds, "sample.id")
refFracPlot(gds, variant.id=5:6, 
            highlight=list(samples[2:3], samples[4:5]))
seqClose(gds)
}

\keyword{manip}
