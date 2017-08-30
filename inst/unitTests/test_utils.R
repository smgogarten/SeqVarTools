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
