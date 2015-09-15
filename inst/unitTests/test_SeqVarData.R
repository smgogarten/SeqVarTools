.testSampleData <- function(gds) {
    require(Biobase)
    sample.id <- seqGetData(gds, "sample.id")
    df <- data.frame(sample.id=seqGetData(gds, "sample.id"),
                     sex=sample(c("M","F"), replace=TRUE, length(sample.id)),
                     stringsAsFactors=FALSE)
    AnnotatedDataFrame(df)
}

test_SeqVarData <- function() {
    ## test creation from existing object
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    sample.id <- seqGetData(gds, "sample.id")
    adf <- .testSampleData(gds)
    svd <- SeqVarData(gds, sampleData=adf)
    checkIdentical(sample.id, seqGetData(svd, "sample.id"))
    seqClose(svd)

    ## test creation from filename
    svd <- SeqVarData(gds.fn, sampleData=adf)
    checkIdentical(sample.id, seqGetData(svd, "sample.id"))
    seqClose(svd)

    ## test errors
    checkException(SeqVarData("abc"))
    checkException(SeqVarData(gds, pData(adf)))
    adf$sample.id <- paste0("a", adf$sample.id)
    checkException(SeqVarData(gds, adf))
    adf$sample.id <- NULL
    checkException(SeqVarData(gds, adf))
}

test_filters <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    adf <- .testSampleData(gds)
    svd <- SeqVarData(gds, sampleData=adf)
    seqSetFilter(svd, sample.id=adf$sample.id[1:10])
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    seqClose(svd)
}

test_replacement <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    adf <- .testSampleData(gds)
    svd <- SeqVarData(gds, sampleData=adf)
    adf$tmp <- "a"
    sampleData(svd) <- adf
    checkIdentical(adf$tmp, sampleData(svd)$tmp)
    checkException(sampleData(svd) <- adf[1:10,])
    seqClose(svd)
}
