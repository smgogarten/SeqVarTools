.testSampleData <- function(gds) {
    require(Biobase)
    sample.id <- seqGetData(gds, "sample.id")
    df <- data.frame(sample.id=sample.id,
                     sex=sample(c("M","F"), replace=TRUE, length(sample.id)),
                     stringsAsFactors=FALSE)
    AnnotatedDataFrame(df)
}

.testVariantData <- function(gds) {
    require(Biobase)
    variant.id <- seqGetData(gds, "variant.id")
    df <- data.frame(variant.id=variant.id,
                     weight=runif(length(variant.id)),
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
    gds <- seqOpen(gds.fn)
    checkException(SeqVarData(gds, pData(adf)))
    adf$sample.id <- paste0("a", adf$sample.id)
    checkException(SeqVarData(gds, adf))
    adf$sample.id <- NULL
    checkException(SeqVarData(gds, adf))
    seqClose(gds)
}

test_variantData <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    adf <- .testVariantData(gds)
    svd <- SeqVarData(gds, variantData=adf)
    checkIdentical(adf$variant.id, seqGetData(svd, "variant.id"))

    ## test errors
    adf$variant.id <- adf$variant.id + 1
    checkException(SeqVarData(gds, variantData=adf))
    
    seqClose(svd)
}

test_filters <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    sdf <- .testSampleData(gds)
    vdf <- .testVariantData(gds)
    svd <- SeqVarData(gds, sampleData=sdf, variantData=vdf)
    seqSetFilter(svd, sample.id=sdf$sample.id[1:10], variant.id=vdf$variant.id[1:10], verbose=FALSE)
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    checkEquals(seqGetData(svd, "variant.id"), variantData(svd)$variant.id)
    seqClose(svd)
}

test_replacement <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    sdf <- .testSampleData(gds)
    vdf <- .testVariantData(gds)
    svd <- SeqVarData(gds, sampleData=sdf, variantData=vdf)
    
    sdf$tmp <- "a"
    sampleData(svd) <- sdf
    checkIdentical(sdf$tmp, sampleData(svd)$tmp)
    checkException(sampleData(svd) <- sdf[1:10,])

    vdf$tmp <- "a"
    variantData(svd) <- vdf
    checkIdentical(vdf$tmp, variantData(svd)$tmp)
    checkException(variantData(svd) <- vdf[1:10,])
    
    seqClose(svd)
}

test_missing <- function() {
    ## test creation with missing arguments
    gds.fn <- seqExampleFileName("gds")
    svd <- SeqVarData(gds.fn)
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    checkEquals(seqGetData(svd, "variant.id"), variantData(svd)$variant.id)
    seqClose(svd)
}

test_granges <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    svd <- SeqVarData(gds)
    checkEquals(granges(gds), granges(svd))
    
    vdf <- .testVariantData(gds)
    variantData(svd) <- vdf
    checkEquals(pData(vdf)[,-1,drop=FALSE], as.data.frame(mcols(granges(svd))))
    
    seqClose(svd)
}

test_filtered <- function() {
    gds.fn <- seqExampleFileName("gds")
    gds <- seqOpen(gds.fn)
    seqSetFilter(gds, sample.sel=1:10, variant.sel=1:10, verbose=FALSE)
    svd <- SeqVarData(gds)
    checkEquals(10, nrow(sampleData(svd)), checkNames=FALSE)
    checkEquals(10, nrow(variantData(svd)), checkNames=FALSE)
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    checkEquals(seqGetData(svd, "variant.id"), variantData(svd)$variant.id)
    seqClose(svd)
}
