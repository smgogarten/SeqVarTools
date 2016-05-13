test_nSamp <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- sum(seqGetFilter(gds)$sample.sel)
  checkEquals(n, SeqVarTools:::.nSamp(gds))

  n <- 10
  seqSetFilter(gds, sample.sel=1:n)
  checkEquals(n, SeqVarTools:::.nSamp(gds))
  seqClose(gds)
}

test_nVar <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- sum(seqGetFilter(gds)$variant.sel)
  checkEquals(n, SeqVarTools:::.nVar(gds))

  n <- 10
  seqSetFilter(gds, variant.sel=1:n)
  checkEquals(n, SeqVarTools:::.nVar(gds))
  seqClose(gds)
}

test_nSampUnfiltered <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- SeqVarTools:::.nSamp(gds)
  seqSetFilter(gds, sample.sel=1:10)
  checkEquals(n, SeqVarTools:::.nSampUnfiltered(gds))
  seqClose(gds)
}

test_nVarUnfiltered <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  n <- SeqVarTools:::.nVar(gds)
  seqSetFilter(gds, sample.sel=1:10)
  checkEquals(n, SeqVarTools:::.nVarUnfiltered(gds))
  seqClose(gds)
}
