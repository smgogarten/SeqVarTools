## constructor
SeqVarData <- function(gds, sampleData, variantData) {
    closeOnFail <- FALSE
    if (is.character(gds) & length(gds) == 1) {
        gds <- seqOpen(gds)
        closeOnFail <- TRUE
    } else if (!is(gds, "SeqVarGDSClass")) {
        stop("gds must be a SeqVarGDSClass object or a valid filename")
    }
    class(gds) <- "SeqVarData"

    if (missing(sampleData)) {
        filt.orig <- seqGetFilter(gds)
        seqSetFilter(gds, sample.sel=1:.nSampUnfiltered(gds), verbose=FALSE)
        sampleData <- AnnotatedDataFrame(data.frame(sample.id=seqGetData(gds, "sample.id"),
                                                    stringsAsFactors=FALSE))
        seqSetFilter(gds, sample.sel=filt.orig$sample.sel, verbose=FALSE)
        
    }
    gds@sampleData <- sampleData
    
    if (missing(variantData)) {
        filt.orig <- seqGetFilter(gds)
        seqSetFilter(gds, variant.sel=1:.nVarUnfiltered(gds), verbose=FALSE)
        variantData <- AnnotatedDataFrame(data.frame(variant.id=seqGetData(gds, "variant.id"),
                                                     stringsAsFactors=FALSE))
        seqSetFilter(gds, variant.sel=filt.orig$variant.sel, verbose=FALSE)
    }
    gds@variantData <- variantData
    
    check <- validObject(gds, test=TRUE)
    if (is.character(check)) {
        if (closeOnFail) seqClose(gds)
        stop(check)
    }
    gds
}

setValidity("SeqVarData",
            function(object) {
                if (!("sample.id" %in% varLabels(sampleData(object)))) {
                    return("sampleData does not have column sample.id")
                }
                if (!identical(seqGetData(object, "sample.id"),
                               sampleData(object)$sample.id)) {
                    return("sample.id in sampleData does not match GDS file")
                }
                if (!("variant.id" %in% varLabels(variantData(object)))) {
                    return("variantData does not have column variant.id")
                }
                if (!identical(seqGetData(object, "variant.id"),
                               variantData(object)$variant.id)) {
                    return("variant.id in variantData does not match GDS file")
                }
                TRUE
            })

            
setMethod("sampleData",
          "SeqVarData",
          function(x) {
              x@sampleData[seqGetFilter(x)$sample.sel,]
          })

setReplaceMethod("sampleData",
                 c("SeqVarData", "AnnotatedDataFrame"),
                 function(x, value) {
                     if (!identical(x@sampleData$sample.id,
                                    value$sample.id)) {
                         stop("sample.id does not match GDS file")
                     }
                     slot(x, "sampleData") <- value
                     x
                 })

            
setMethod("variantData",
          "SeqVarData",
          function(x) {
              x@variantData[seqGetFilter(x)$variant.sel,]
          })

setReplaceMethod("variantData",
                 c("SeqVarData", "AnnotatedDataFrame"),
                 function(x, value) {
                     if (!identical(x@variantData$variant.id,
                                    value$variant.id)) {
                         stop("variant.id does not match GDS file")
                     }
                     slot(x, "variantData") <- value
                     x
                 })


setMethod("granges",
          "SeqVarData",
          function(x, ...) {
              gr <- callNextMethod(x, ...)
              df <- pData(variantData(x))
              df$variant.id <- NULL # redundant with names of gr
              mcols(gr) <- cbind(mcols(gr), as(df, "DataFrame"))
              gr
          })
