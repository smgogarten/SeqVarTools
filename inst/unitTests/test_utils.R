test_rangesTo <- function() {
    gds <- SeqVarTools:::.testData()
    ranges <- GRanges("1", IRanges(start=c(1e6,3e6), end=c(2e6,4e6)))
    checkEquals(1:7, SeqVarTools:::.rangesToSel(gds, ranges))
    checkEquals(1:7, SeqVarTools:::.rangesToID(gds, ranges))

    # check empty ranges
    ranges <- GRanges("X", IRanges(start=c(1e6,3e6), end=c(2e6,4e6)))
    checkEquals(0, length(SeqVarTools:::.rangesToSel(gds, ranges)))
    seqClose(gds)
}

test_nSamp <- function() {
  gds <- SeqVarTools:::.testData()
  n <- sum(seqGetFilter(gds)$sample.sel)
  checkEquals(n, SeqVarTools:::.nSamp(gds))

  n <- 10
  seqSetFilter(gds, sample.sel=1:n)
  checkEquals(n, SeqVarTools:::.nSamp(gds))
  seqClose(gds)
}

test_nVar <- function() {
  gds <- SeqVarTools:::.testData()
  n <- sum(seqGetFilter(gds)$variant.sel)
  checkEquals(n, SeqVarTools:::.nVar(gds))

  n <- 10
  seqSetFilter(gds, variant.sel=1:n)
  checkEquals(n, SeqVarTools:::.nVar(gds))
  seqClose(gds)
}

test_nSampUnfiltered <- function() {
  gds <- SeqVarTools:::.testData()
  n <- SeqVarTools:::.nSamp(gds)
  seqSetFilter(gds, sample.sel=1:10)
  checkEquals(n, SeqVarTools:::.nSampUnfiltered(gds))
  seqClose(gds)
}

test_nVarUnfiltered <- function() {
  gds <- SeqVarTools:::.testData()
  n <- SeqVarTools:::.nVar(gds)
  seqSetFilter(gds, sample.sel=1:10)
  checkEquals(n, SeqVarTools:::.nVarUnfiltered(gds))
  seqClose(gds)
}

test_emptyGenoMatrix <- function() {
    gds <- SeqVarTools:::.testData()
    SeqVarTools:::.emptyVarFilter(gds)
    checkTrue(SeqVarTools:::.emptyDim(gds))
    m <- SeqVarTools:::.emptyGenoMatrix(gds, use.names=TRUE)
    checkEquals(SeqVarTools:::.nSamp(gds), nrow(m))
    checkEquals(0, ncol(m))
    checkEquals(seqGetData(gds, "sample.id"), rownames(m))

    seqResetFilter(gds, verbose=FALSE)
    SeqVarTools:::.emptySampFilter(gds)
    checkTrue(SeqVarTools:::.emptyDim(gds))
    m <- SeqVarTools:::.emptyGenoMatrix(gds, use.names=TRUE)
    checkEquals(SeqVarTools:::.nVar(gds), ncol(m))
    checkEquals(0, nrow(m))
    checkEquals(as.character(seqGetData(gds, "variant.id")), colnames(m))

    seqClose(gds)
}

test_nSampObserved <- function() {
    gds <- SeqVarTools:::.testData()
    nso <- SeqVarTools:::.nSampObserved(gds)
    checkEquals(SeqVarTools:::.nVar(gds), length(nso))
    checkEquals(SeqVarTools:::.nSamp(gds), max(nso))
    seqClose(gds)
}

test_nSamp_empty <- function() {
    gds <- SeqVarTools:::.testData()
    seqSetFilter(gds, sample.sel=rep(SeqVarTools:::.nSamp(gds), FALSE), verbose=FALSE)
    checkEquals(0, SeqVarTools:::.nSamp(gds))
    checkEquals(rep(0, SeqVarTools:::.nVar(gds)), SeqVarTools:::.nSampObserved(gds))
    seqClose(gds)
}
