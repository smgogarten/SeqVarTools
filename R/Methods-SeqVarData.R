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
        sampleData <- AnnotatedDataFrame(data.frame(matrix(nrow=.nSampUnfiltered(gds), ncol=0)))
    }
    gds@sampleData <- sampleData
    
    if (missing(variantData)) {
        variantData <- AnnotatedDataFrame(data.frame(matrix(nrow=.nVarUnfiltered(gds), ncol=0)))
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
                if (ncol(sampleData(object)) > 0) {
                    if (!("sample.id" %in% varLabels(sampleData(object)))) {
                        return("sampleData does not have column sample.id")
                    }
                    if (!identical(seqGetData(object, "sample.id"),
                                   sampleData(object)$sample.id)) {
                        return("sample.id in sampleData does not match GDS file")
                    }
                }
                if (ncol(variantData(object)) > 0) {
                    if (!("variant.id" %in% varLabels(variantData(object)))) {
                        return("variantData does not have column variant.id")
                    }
                    if (!identical(seqGetData(object, "variant.id"),
                                   variantData(object)$variant.id)) {
                        return("variant.id in variantData does not match GDS file")
                    }
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
                     seqSetFilter(x, sample.sel=1:.nSampUnfiltered(x), action="push+set", verbose=FALSE)
                     if (!identical(seqGetData(x, "sample.id"), value$sample.id)) {
                         stop("sample.id does not match GDS file")
                     }
                     seqSetFilter(x, action="pop", verbose=FALSE)
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
                     seqSetFilter(x, variant.sel=1:.nVarUnfiltered(x), action="push+set", verbose=FALSE)
                     if (!identical(seqGetData(x, "variant.id"), value$variant.id)) {
                         stop("variant.id does not match GDS file")
                     }
                     seqSetFilter(x, action="pop", verbose=FALSE)
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
