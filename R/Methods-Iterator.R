SeqVarIterator <- function(seqData, variantFilter, verbose=TRUE) {
    class(seqData) <- "SeqVarIterator"
    seqData@variantFilter <- variantFilter

    ## set filter to first element
    seqSetFilter(seqData, variant.sel=variantFilter[[1]], verbose=verbose)
    
    ## pass-by-reference slot for lastFilter
    seqData@lastFilter <- new.env()
    lastFilter(seqData) <- 1
    
    seqData
}


SeqVarBlockIterator <- function(seqData, variantBlock=10000, verbose=TRUE) {

    if (variantBlock > .nVar(seqData)) {
        variantBlock <- .nVar(seqData)
    }
    
    ## original filter
    current.sel <- which(seqGetFilter(seqData)$variant.sel)

    ## divide selected variants into blocks
    variantFilter <- split(current.sel, ceiling(seq_along(current.sel)/variantBlock))
    
    seqData <- SeqVarIterator(seqData, variantFilter, verbose=verbose)
    class(seqData) <- "SeqVarBlockIterator"
    seqData@variantBlock <- as.integer(variantBlock)

    seqData
}


SeqVarRangeIterator <- function(seqData, variantRanges, verbose=TRUE) {
    ## original filter
    current.sel <- which(seqGetFilter(seqData)$variant.sel)

    ## filter for iterator ranges
    variants <- granges(seqData)
    new.sel <- .subjectByQuery(variantRanges, variants, hits.only=FALSE)$subjectHits
    variantFilter <- lapply(new.sel, function(x) current.sel[x])
    
    seqData <- SeqVarIterator(seqData, variantFilter, verbose=verbose)
    class(seqData) <- "SeqVarRangeIterator"
    seqData@variantRanges <- variantRanges

    seqData
}


SeqVarWindowIterator <- function(seqData, windowSize=10000, windowShift=5000, verbose=TRUE) {
    ## original filter
    current.sel <- which(seqGetFilter(seqData)$variant.sel)
    
    # identify windows
    variants <- granges(seqData)
    windows <- do.call("c", lapply(unique(seqnames(variants)), function(chr) {
        start <- seq(1, max(end(variants[seqnames(variants) == chr])), windowShift)
        GRanges(seqnames=chr, ranges=IRanges(start=start, width=windowSize))
    }))

    ## variants in each window
    sbq <- .subjectByQuery(windows, variants, hits.only=TRUE)
    new.sel <- sbq$subjectHits

    ## only keep unique windows
    keep <- !duplicated(new.sel)
    new.sel <- new.sel[keep]
    windows <- windows[sbq$queryHits[keep]]

    variantFilter <- lapply(new.sel, function(x) current.sel[x])
    
    seqData <- SeqVarIterator(seqData, variantFilter, verbose=verbose)
    class(seqData) <- "SeqVarWindowIterator"
    seqData@variantRanges <- windows
    seqData@windowSize <- as.integer(windowSize)
    seqData@windowShift <- as.integer(windowShift)
        
    seqData
}


SeqVarListIterator <- function(seqData, variantRanges, verbose=TRUE) {
    ## original filter
    current.sel <- which(seqGetFilter(seqData)$variant.sel)
    
    ## filter for iterator ranges
    variants <- granges(seqData)
    new.sel <- .subjectByQuery(variantRanges, variants, hits.only=FALSE)$subjectHits
    variantFilter <- lapply(new.sel, function(x) current.sel[x])
    
    seqData <- SeqVarIterator(seqData, variantFilter, verbose=verbose)
    class(seqData) <- "SeqVarListIterator"
    seqData@variantRanges <- variantRanges

    seqData
}


setMethod("variantFilter",
          "SeqVarIterator",
          function(x) {
              x@variantFilter
          })

setMethod("lastFilter",
          "SeqVarIterator",
          function(x) {
              x@lastFilter$i
          })

setReplaceMethod("lastFilter",
                 c("SeqVarIterator", "numeric"),
                 function(x, value) {
                     x@lastFilter$i <- as.integer(value)
                     x
                 })

setMethod("iterateFilter",
          "SeqVarIterator",
          function(x, verbose=TRUE) {
              ## set filter for next element
              if (lastFilter(x) < length(variantFilter(x))) {
                  i <- lastFilter(x) + 1
                  seqSetFilter(x, variantFilter(x)[[i]], verbose=verbose)
                  lastFilter(x) <- i
                  return(TRUE)
              } else {
                  .emptyVarFilter(x, verbose=verbose)
                  return(FALSE)
              }
          })


setMethod("variantBlock",
          "SeqVarBlockIterator",
          function(x) {
              x@variantBlock
          })


setMethod("variantRanges",
          "SeqVarRangeIterator",
          function(x) {
              x@variantRanges
          })

setMethod("variantRanges",
          "SeqVarListIterator",
          function(x) {
              x@variantRanges
          })
