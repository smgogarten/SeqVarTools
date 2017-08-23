
test_countSingletons <- function() {
  gds <- SeqVarTools:::.testData()
  geno <- altDosage(gds)
  hasalt <- geno > 0 & !is.na(geno)
  singleton <- colSums(hasalt) == 1
  singmat <- matrix(singleton, nrow=nrow(geno), ncol=ncol(geno), byrow=TRUE)
  ct <- rowSums(hasalt & singmat)
  checkEquals(ct, countSingletons(gds), checkNames=FALSE)
  seqClose(gds)
}
