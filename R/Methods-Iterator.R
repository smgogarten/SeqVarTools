
SeqVarBlockIterator <- function(seqData, variantBlock=10000, verbose=TRUE) {
    class(seqData) <- "SeqVarBlockIterator"
    seqData@variantBlock <- as.integer(variantBlock)

    ## set filter to first block (store original with push to stack)
    seqSetFilter(seqData, variant.sel=1:variantBlock, action="push+intersect", verbose=verbose)

    ## pass-by-reference slot for lastVariant
    seqData@lastVariant <- new.env()
    lastVariant(seqData) <- variantBlock

    seqData
}


SeqVarRangeIterator <- function(seqData, variantRanges, verbose=TRUE) {
    class(seqData) <- "SeqVarRangeIterator"
    seqData@variantRanges <- variantRanges

    ## store original filter
    seqSetFilter(seqData, action="push", verbose=verbose)
    
    ## set filter to first range
    seqSetFilter(seqData, variantRanges[1], intersect=TRUE, verbose=verbose)

    ## pass-by-reference slot for lastRange
    seqData@lastRange <- new.env()
    lastRange(seqData) <- 1

    seqData
}


SeqVarWindowIterator <- function(seqData, windowSize=10000, windowShift=5000, verbose=TRUE) {
    # identify windows
    variants <- granges(seqData)
    windows <- do.call("c", lapply(unique(seqnames(variants)), function(chr) {
        start <- seq(1, max(end(variants[seqnames(variants) == chr])), windowShift)
        windows.chr <- GRanges(seqnames=chr, ranges=IRanges(start=start, width=windowSize))
        .subsetByUniqueOverlaps(query=windows.chr, subject=variants)
    }))

    seqData <- SeqVarRangeIterator(seqData, variantRanges=windows, verbose=verbose)
    class(seqData) <- "SeqVarWindowIterator"
    seqData@windowSize <- as.integer(windowSize)
    seqData@windowShift <- as.integer(windowShift)
        
    seqData
}


setMethod("restoreFilter",
          "SeqVarGDSClass",
          function(x, verbose=TRUE) {
              seqSetFilter(x, action="pop", verbose=verbose)
          })

setMethod("variantBlock",
          "SeqVarBlockIterator",
          function(x) {
              x@variantBlock
          })

setMethod("lastVariant",
          "SeqVarBlockIterator",
          function(x) {
              x@lastVariant$i
          })

setReplaceMethod("lastVariant",
                 c("SeqVarBlockIterator", "numeric"),
                 function(x, value) {
                     x@lastVariant$i <- as.integer(value)
                     x
                 })

setMethod("iterateFilter",
          "SeqVarBlockIterator",
          function(x, verbose=TRUE) {
              ## restore original filter
              restoreFilter(x, verbose=verbose)
              nVar <- .nVar(x)

              ## set filter for next block
              if (lastVariant(x) < nVar) {
                  start <- lastVariant(x) + 1
                  end <- min(lastVariant(x) + variantBlock(x), nVar)
                  seqSetFilter(x, variant.sel=start:end, action="push+intersect", verbose=verbose)
                  lastVariant(x) <- end
                  return(TRUE)
              } else {
                  .emptyVarFilter(x, verbose=verbose)
                  return(FALSE)
              }
          })


setMethod("variantRanges",
          "SeqVarRangeIterator",
          function(x) {
              x@variantRanges
          })

setMethod("lastRange",
          "SeqVarRangeIterator",
          function(x) {
              x@lastRange$i
          })

setReplaceMethod("lastRange",
                 c("SeqVarRangeIterator", "numeric"),
                 function(x, value) {
                     x@lastRange$i <- as.integer(value)
                     x
                 })

setMethod("iterateFilter",
          "SeqVarRangeIterator",
          function(x, verbose=TRUE) {
              ## restore original filter
              restoreFilter(x, verbose=verbose)

              ## set filter for next range
              if (lastRange(x) < length(variantRanges(x))) {
                  i <- lastRange(x) + 1
                  seqSetFilter(x, action="push", verbose=verbose)
                  seqSetFilter(x, variantRanges(x)[i], intersect=TRUE, verbose=verbose)
                  lastRange(x) <- i
                  return(TRUE)
              } else {
                  .emptyVarFilter(x, verbose=verbose)
                  return(FALSE)
              }
          })

