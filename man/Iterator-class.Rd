\name{Iterator}
\docType{class}

\alias{SeqVarIterator}
\alias{SeqVarIterator-class}
\alias{SeqVarBlockIterator}
\alias{SeqVarBlockIterator-class}
\alias{SeqVarRangeIterator}
\alias{SeqVarRangeIterator-class}
\alias{SeqVarWindowIterator}
\alias{SeqVarWindowIterator-class}
\alias{SeqVarListIterator}
\alias{SeqVarListIterator-class}
\alias{show,SeqVarIterator-method}
\alias{show,SeqVarRangeIterator-method}
\alias{show,SeqVarListIterator-method}
\alias{iterateFilter}
\alias{iterateFilter,SeqVarIterator-method}
\alias{lastFilter}
\alias{lastFilter,SeqVarIterator-method}
\alias{lastFilter<-}
\alias{lastFilter<-,SeqVarIterator,numeric-method}
\alias{variantFilter}
\alias{variantFilter,SeqVarIterator-method}
\alias{variantBlock}
\alias{variantBlock,SeqVarBlockIterator-method}
\alias{variantRanges}
\alias{variantRanges,SeqVarRangeIterator-method}
\alias{variantRanges,SeqVarListIterator-method}
\alias{currentRanges}
\alias{currentRanges,SeqVarBlockIterator-method}
\alias{currentRanges,SeqVarRangeIterator-method}
\alias{currentRanges,SeqVarListIterator-method}

\title{Iterators}

\description{
	Extends \code{\link{SeqVarData}} to provide iterators over variants.
}

\details{
	Iterator classes allow for iterating filters over blocks of variants, ranges, or sliding windows.

	For \code{SeqVarBlockIterator}, each call to \code{iterateFilter} will select the next unit of \code{variantBlock} variants.

	For \code{SeqVarRangeIterator}, each call to \code{iterateFilter} will select the next range in \code{variantRanges}.

	\code{SeqVarWindowIterator} is an extension of \code{SeqVarRangeIterator} where the ranges are determined automatically by sliding a window of size \code{windowSize} base pairs by steps of \code{windowShift} across the genome. Only windows containing unique sets of variants are kept.

	For \code{SeqVarListIterator}, each call to \code{iterateFilter} will select the next set of ranges in \code{variantRanges}.

	Any filter set on the object previously will be applied in addition to the selected blocks or ranges. The original filter can be restored after iteration by calling \code{restoreFilter}.
}


\section{Constructors}{
  \describe{
    \item{}{
      \code{SeqVarBlockIterator(seqData, variantBlock=10000, verbose=TRUE)}: Returns a \code{SeqVarBlockIterator} object with the filter set to the first block.

      \code{seqData} is a \code{\link{SeqVarData}} object.

      \code{variantBlock} is an integer specifying the number of variants in an iteration block.

      \code{verbose} is a logical indicator for verbose output.
    }
    \item{}{
      \code{SeqVarRangeIterator(seqData, variantRanges=GRanges(), verbose=TRUE)}: Returns a \code{SeqVarRangeIterator} object with the filter set to the first range.

      \code{seqData} is a \code{\link{SeqVarData}} object.

      \code{variantRanges} is a \code{\link{GRanges}} object specifying the ranges for iteration.

      \code{verbose} is a logical indicator for verbose output.
    }
    \item{}{
      \code{SeqVarWindowIterator(seqData,  windowSize=10000, windowShift=5000, verbose=TRUE)}: Returns a \code{SeqVarWindowIterator} object with the filter set to the first window.

      \code{seqData} is a \code{\link{SeqVarData}} object.

      \code{windowSize} is the size in base pairs of the sliding window.

      \code{windowShift} is the size in base pairs of the step for each consecutive window.

      \code{verbose} is a logical indicator for verbose output.
    }
    \item{}{
      \code{SeqVarListIterator(seqData, variantRanges, verbose=TRUE)}: Returns a \code{SeqVarRangeIterator} object with the filter set to the first range.

      \code{seqData} is a \code{\link{SeqVarData}} object.

      \code{variantRanges} is a \code{\link{GRangesList}} object specifying the ranges for iteration.

      \code{verbose} is a logical indicator for verbose output.
    }
  }
}


\section{Accessors}{
  \describe{
    \item{}{\code{iterateFilter(x)}: Advance the filter to the next block, range, or set of ranges. Returns \code{TRUE} while there are still variants left to be read, \code{FALSE} if the end of iteration is reached.
    }
    \item{}{\code{lastFilter(x)}, \code{lastFilter(x)<- value}: Get or set the last filter index from the previous call to \code{iterateFilter}.
    }
    \item{}{\code{variantBlock(x)}: Get the size of the variant block.
    }
    \item{}{\code{variantFilter(x)}: Get the list of variant indices.
    }
    \item{}{\code{variantRanges(x)}: Get the variant ranges.
    }
    \item{}{\code{currentRanges(x)}: Get the ranges selected in the current iteration.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{SeqVarData}},
  \code{\link{seqSetFilter}}
}

\examples{
gds <- seqOpen(seqExampleFileName("gds"))
seqData <- SeqVarData(gds)

# iterate by blocks
seqSetFilter(seqData, variant.sel=seq(1,1000,2))
iterator <- SeqVarBlockIterator(seqData, variantBlock=10)
seqGetData(iterator, "variant.id")
iterateFilter(iterator)
seqGetData(iterator, "variant.id")
seqResetFilter(iterator)

# iterate by ranges
library(GenomicRanges)
gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
iterator <- SeqVarRangeIterator(seqData, variantRanges=gr)
granges(iterator)
iterateFilter(iterator) # no variants in the second range
granges(iterator)
iterateFilter(iterator)
granges(iterator)
iterateFilter(iterator)
seqResetFilter(iterator)

# iterate by windows
seqSetFilterChrom(seqData, include="22")
iterator <- SeqVarWindowIterator(seqData)
seqGetData(iterator, "variant.id")
while (iterateFilter(iterator)) {
    print(seqGetData(iterator, "variant.id"))
}
seqResetFilter(iterator)

# iterate by list of ranges
gr <- GRangesList(
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(16e6, 17e6), width=1e6)),
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(18e6, 20e6), width=1e6)))
iterator <- SeqVarListIterator(seqData, variantRanges=gr)
granges(iterator)
iterateFilter(iterator)
granges(iterator)
iterateFilter(iterator)
seqResetFilter(iterator)

seqClose(iterator)
}

\keyword{methods}
\keyword{classes}
