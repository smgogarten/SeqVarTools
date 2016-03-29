test_nSamp <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- sum(seqGetFilter(gds)$sample.sel)
  checkEquals(n, SeqVarTools:::.nSamp(gds))

  n <- 10
  seqSetFilter(gds, sample.id=seqGetData(gds, "sample.id")[1:n])
  checkEquals(10, SeqVarTools:::.nSamp(gds))
  seqClose(gds)
}

test_nVar <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- sum(seqGetFilter(gds)$variant.sel)
  checkEquals(n, SeqVarTools:::.nVar(gds))

  n <- 10
  seqSetFilter(gds, variant.id=1:n)
  checkEquals(10, SeqVarTools:::.nVar(gds))
  seqClose(gds)
}
