test_pca <- function() {
  gds <- SeqVarTools:::.testData()
  pc <- pca(gds, eigen.cnt=8)
  checkEquals(8, length(pc$eigenval))
  checkEquals(c(length(seqGetData(gds, "sample.id")), 8),
              dim(pc$eigenvect))
  seqClose(gds)
}

test_pca_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  pc <- pca(gds)
  seqSetFilter(gds)
  checkIdentical(pc,
                 applyMethod(gds, pca, variant=var.id, sample=samp.id))
  seqClose(gds)
}

