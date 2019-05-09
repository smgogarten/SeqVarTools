.testSampleData <- function(gds) {
    require(Biobase)
    sample.id <- seqGetData(gds, "sample.id")
    set.seed(66); sex <- sample(c("M","F"), replace=TRUE, length(sample.id))
    df <- data.frame(sample.id, sex, stringsAsFactors=FALSE)
    AnnotatedDataFrame(df)
}

.testVariantData <- function(gds) {
    require(Biobase)
    variant.id <- seqGetData(gds, "variant.id")
    set.seed(77); weight <- runif(length(variant.id))
    df <- data.frame(variant.id, weight, stringsAsFactors=FALSE)
    AnnotatedDataFrame(df)
}

test_SeqVarData <- function() {
    ## test creation from existing object
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
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
    gds <- SeqVarTools:::.testData()
    adf <- .testVariantData(gds)
    svd <- SeqVarData(gds, variantData=adf)
    checkIdentical(adf$variant.id, seqGetData(svd, "variant.id"))

    ## test errors
    adf$variant.id <- adf$variant.id + 1
    checkException(SeqVarData(gds, variantData=adf))
    
    seqClose(svd)
}

test_filters <- function() {
    gds <- SeqVarTools:::.testData()
    sdf <- .testSampleData(gds)
    vdf <- .testVariantData(gds)
    svd <- SeqVarData(gds, sampleData=sdf, variantData=vdf)
    seqSetFilter(svd, sample.id=sdf$sample.id[1:10], variant.id=vdf$variant.id[1:10], verbose=FALSE)
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    checkEquals(seqGetData(svd, "variant.id"), variantData(svd)$variant.id)
    seqClose(svd)
}

test_replacement <- function() {
    gds <- SeqVarTools:::.testData()
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
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gds.fn <- seqExampleFileName("gds")
    svd <- SeqVarData(gds.fn)
    checkEquals(length(seqGetData(svd, "sample.id")), nrow(sampleData(svd)), checkNames=FALSE)
    checkEquals(length(seqGetData(svd, "variant.id")), nrow(variantData(svd)), checkNames=FALSE)
    seqClose(svd)
}

test_granges <- function() {
    gds <- SeqVarTools:::.testData()
    svd <- SeqVarData(gds)
    checkEquals(granges(gds), granges(svd))
    
    vdf <- .testVariantData(gds)
    variantData(svd) <- vdf
    x <- pData(vdf)[,-1,drop=FALSE]
    row.names(x) <- as.character(row.names(x))
    checkEquals(x, as.data.frame(S4Vectors::mcols(granges(svd))), checkNames=FALSE)
    
    seqClose(svd)
}

test_filtered <- function() {
    gds <- SeqVarTools:::.testData()
    sdf <- .testSampleData(gds)
    vdf <- .testVariantData(gds)
    seqSetFilter(gds, sample.sel=1:10, variant.sel=1:10, verbose=FALSE)
    svd <- SeqVarData(gds)
    checkEquals(10, nrow(sampleData(svd)), checkNames=FALSE)
    checkEquals(10, nrow(variantData(svd)), checkNames=FALSE)
    svd <- SeqVarData(gds, sampleData=sdf, variantData=vdf)
    checkEquals(seqGetData(svd, "sample.id"), sampleData(svd)$sample.id)
    checkEquals(seqGetData(svd, "variant.id"), variantData(svd)$variant.id)
    seqClose(svd)
}

test_validateSex <- function() {
    svd <- SeqVarTools:::.testSeqVarData()
    checkTrue(is.null(validateSex(svd)))
    sdf <- .testSampleData(svd)
    sampleData(svd) <- sdf
    checkTrue(!is.null(validateSex(svd)))
    set.seed(88); sampleData(svd)$sex <- sample(c(1,2,NA), nrow(sdf), replace=TRUE)
    checkTrue(setequal(c("M","F",NA), validateSex(svd)))
    seqClose(svd)
}
