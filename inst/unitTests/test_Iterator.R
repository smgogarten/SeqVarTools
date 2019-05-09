library(GenomicRanges)

test_iterator_block <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    it <- SeqVarBlockIterator(gds, variantBlock=100, verbose=FALSE)
    checkEquals(1:100, seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(101:200, seqGetData(it, "variant.id"))
    checkEquals(granges(it), currentRanges(it))
    seqClose(it)
}

test_iterator_block_large <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    it <- SeqVarBlockIterator(gds, variantBlock=10000, verbose=FALSE)
    checkEquals(1:SeqVarTools:::.nVar(it), seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(GRanges(), currentRanges(it), checkNames=FALSE)
    seqClose(it)
}

test_iterator_block_prev <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    set.seed(11); var.sel <- sort(sample(1:SeqVarTools:::.nVar(gds), 250))
    seqSetFilter(gds, variant.sel=var.sel, verbose=FALSE)
    it <- SeqVarBlockIterator(gds, variantBlock=100, verbose=FALSE)
    checkEquals(var.sel[1:100], seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(var.sel[101:200], seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(var.sel[201:250], seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_block_samples <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    set.seed(22); samp.sel <- sort(sample(1:SeqVarTools:::.nSamp(gds), 50))
    seqSetFilter(gds, sample.sel=samp.sel, verbose=FALSE)
    it <- SeqVarBlockIterator(gds, variantBlock=100, verbose=FALSE)
    checkEquals(samp.sel, which(seqGetFilter(it)$sample.sel))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(samp.sel, which(seqGetFilter(it)$sample.sel))
    seqClose(it)
}

test_iterator_range <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    it <- SeqVarRangeIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(1:3, seqGetData(it, "variant.id"))
    checkEquals(gr[1], currentRanges(it))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(4:7, seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_range_prev <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    seqSetFilter(gds, variant.sel=c(1,3,5,7), verbose=FALSE)
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    it <- SeqVarRangeIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(c(1,3), seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(c(5,7), seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_range_samples <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    set.seed(33); samp.sel <- sort(sample(1:SeqVarTools:::.nSamp(gds), 50))
    seqSetFilter(gds, sample.sel=samp.sel, verbose=FALSE)
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    it <- SeqVarRangeIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(samp.sel, which(seqGetFilter(it)$sample.sel))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(samp.sel, which(seqGetFilter(it)$sample.sel))
    seqClose(it)
}

test_unique_overlaps <- function() {
    subject <- GRanges(seqnames="1", ranges=IRanges(seq(10,50,10), (seq(10,50,10))))
    query <- GRanges(seqnames="1", ranges=IRanges(c(1,5,9,25,25,29), c(5,15,18,45,55,41)))
    chk <- list(queryHits=2:6, subjectHits=list(1, 1, 3:4, 3:5, 3:4))
    checkEquals(chk, SeqVarTools:::.subjectByQuery(query, subject, hits.only=TRUE))
    chk <- list(queryHits=1:6, subjectHits=list(integer(0), 1, 1, 3:4, 3:5, 3:4))
    checkEquals(chk, SeqVarTools:::.subjectByQuery(query, subject, hits.only=FALSE))
    keep <- c(2,4,5)
    checkEquals(query[keep], SeqVarTools:::.subsetByUniqueOverlaps(query, subject))
}

test_iterator_window <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    seqSetFilterChrom(gds, include="22", verbose=FALSE)
    it <- SeqVarWindowIterator(gds, verbose=FALSE)
    checkTrue(all(width(it@variantRanges) == 10000))
    var <- list(seqGetData(it, "variant.id"))
    i <- 2
    while (iterateFilter(it, verbose=FALSE)) {
        var[[i]] <- seqGetData(it, "variant.id")
        i <- i + 1
    }
    checkEquals(length(var), length(it@variantRanges))
    seqSetFilterChrom(gds, include="22", verbose=FALSE)
    checkEquals(sort(unique(unlist(var))), seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_list <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    gr <- GRangesList(
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6)),
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(3e6, 34e6), width=1e6)))
    it <- SeqVarListIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(1:7, seqGetData(it, "variant.id"))
    checkEquals(gr[[1]], currentRanges(it))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(4:11, seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_list_duplicates <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    gr <- granges(gds)
    grl <- GRangesList(gr[c(1,1,2)])
    it <- SeqVarListIterator(gds, variantRanges=grl, verbose=FALSE)
    checkEquals(1:2, seqGetData(it, "variant.id"))
    seqClose(it)
}


.testVarData <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    vdf <- Biobase::AnnotatedDataFrame(data.frame(
        variant.id=seqGetData(gds, "variant.id"),
        weight=1, stringsAsFactors=FALSE))
    variantData(gds) <- vdf
    gds
}

test_currentRanges_block <- function() {
    gds <- .testVarData()
    it <- SeqVarBlockIterator(gds, variantBlock=10, verbose=FALSE)
    checkEquals(currentRanges(it)$weight, rep(1,10))
    seqClose(it)
}

test_currentRanges_list <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    gr <- GRangesList(
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6), weight=rep(1,2)),
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(3e6, 34e6), width=1e6), weight=rep(2,2)))
    it <- SeqVarListIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(currentRanges(it)$weight, c(1,1))
    iterateFilter(it, verbose=FALSE)
    checkEquals(currentRanges(it)$weight, c(2,2))
    seqClose(it)
}

test_currentRanges_range <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6), weight=rep(1,3))
    it <- SeqVarRangeIterator(gds, variantRanges=gr, verbose=FALSE)
    checkEquals(currentRanges(it)$weight, 1)
    iterateFilter(it, verbose=FALSE)
    checkEquals(currentRanges(it)$weight, 1)
    seqClose(it)
}

test_currentVariants_window <- function() {
    gds <- .testVarData()
    it <- SeqVarWindowIterator(gds, verbose=FALSE)
    checkTrue(all(currentVariants(it)$weight == 1))
    seqClose(it)
}

test_currentVariants_range <- function() {
    gds <- .testVarData()
    gr <- GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6), weight2=rep(2,2))
    it <- SeqVarRangeIterator(gds, variantRanges=gr, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(rownames(cv), as.character(1:3))
    checkEquals(cv$weight, rep(1,3))
    checkEquals(cv$weight2, rep(2,3))
    iterateFilter(it, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(rownames(cv), as.character(4:7))
    checkEquals(cv$weight, rep(1,4))
    checkEquals(cv$weight2, rep(2,4))
    seqClose(it)
}

test_currentVariants_list <- function() {
    gds <- .testVarData()
    gr <- GRangesList(
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6), weight2=rep(2,2)),
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(3e6, 34e6), width=1e6), weight2=rep(3,2)))
    it <- SeqVarListIterator(gds, variantRanges=gr, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(rownames(cv), as.character(1:7))
    checkEquals(cv$weight, rep(1,7))
    checkEquals(cv$weight2, rep(2,7))
    iterateFilter(it, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(rownames(cv), as.character(4:11))
    checkEquals(cv$weight, rep(1,8))
    checkEquals(cv$weight2, rep(3,8))
    seqClose(it)
}

test_currentVariants_block <- function() {
    gds <- .testVarData()
    it <- SeqVarBlockIterator(gds, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(names(cv), c("variant", "weight"))
    checkEquals(rownames(cv), as.character(seqGetData(it, "variant.id")))
    checkTrue(all(cv$weight == 1))
    seqClose(it)
}

test_currentVariants_no_variantData <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    it <- SeqVarWindowIterator(gds, verbose=FALSE)
    cv <- currentVariants(it)
    checkEquals(names(cv), c("variant", "range"))
    seqClose(it)
}

test_resetIterator <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    it <- SeqVarBlockIterator(gds, variantBlock=1000, verbose=FALSE)
    while(iterateFilter(it, verbose=FALSE)) {}
    checkEquals(SeqVarTools:::.nVar(it), 0)
    resetIterator(it, verbose=FALSE)
    checkEquals(SeqVarTools:::.nVar(it), 1000)
    seqClose(it)
}
