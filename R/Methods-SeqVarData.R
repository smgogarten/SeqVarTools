## constructor
SeqVarData <- function(gds, sampleData) {
    closeOnFail <- FALSE
    if (is.character(gds) & length(gds) == 1) {
        gds <- seqOpen(gds)
        closeOnFail <- TRUE
    } else if (!is(gds, "SeqVarGDSClass")) {
        stop("gds must be a SeqVarGDSClass object or a valid filename")
    }
    class(gds) <- "SeqVarData"
    gds@sampleData <- sampleData
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

