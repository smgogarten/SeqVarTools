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
