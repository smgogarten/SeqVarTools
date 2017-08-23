library(GenomicRanges)

test_restore_filter <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    seqSetFilter(gds, variant.sel=1:10, verbose=FALSE)
    checkEquals(1:10, seqGetData(gds, "variant.id"))
    SeqVarTools:::.emptyVarFilter(gds, verbose=FALSE)
    checkEquals(integer(), seqGetData(gds, "variant.id"))
    restoreFilter(gds, verbose=FALSE)
    checkEquals(1:10, seqGetData(gds, "variant.id"))
    SeqVarTools:::.emptyVarFilter(gds, verbose=FALSE)
    restoreFilter(gds, verbose=FALSE)
    checkEquals(1:10, seqGetData(gds, "variant.id"))
    seqClose(gds)
}

test_iterator_block <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    it <- SeqVarBlockIterator(gds, variantBlock=100, verbose=FALSE)
    checkEquals(1:100, seqGetData(it, "variant.id"))
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(101:200, seqGetData(it, "variant.id"))
    seqClose(it)
}

test_iterator_block_prev <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    var.sel <- sort(sample(1:SeqVarTools:::.nVar(gds), 250))
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
    samp.sel <- sort(sample(1:SeqVarTools:::.nSamp(gds), 50))
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
    samp.sel <- sort(sample(1:SeqVarTools:::.nSamp(gds), 50))
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
    restoreFilter(it, verbose=FALSE)
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
    checkTrue(iterateFilter(it, verbose=FALSE))
    checkEquals(4:11, seqGetData(it, "variant.id"))
    checkTrue(!iterateFilter(it, verbose=FALSE))
    checkEquals(integer(), seqGetData(it, "variant.id"))
    seqClose(it)
}
