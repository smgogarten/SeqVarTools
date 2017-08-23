test_alleleFrequency_sum <- function() {
  gds <- SeqVarTools:::.testData()
  maxn <- max(nAlleles(gds))
  af <- matrix(nrow=SeqVarTools:::.nVar(gds), ncol=maxn)
  for (n in 1:maxn) af[,n] <- alleleFrequency(gds, n=(n-1))
  checkTrue(all(rowSums(af) == 1))
  seqClose(gds)
}

test_alleleFrequency_info <- function() {
  gds <- SeqVarTools:::.testData()
  ac <- seqGetData(gds, "annotation/info/AC")
  an <- seqGetData(gds, "annotation/info/AN")
  checkEquals(ac/an, alleleFrequency(gds, n=1))
  seqClose(gds)
}

test_alleleFrequency_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  af <- alleleFrequency(gds)
  seqSetFilter(gds)
  checkIdentical(af,
                 applyMethod(gds, alleleFrequency, variant=var.id, sample=samp.id))
  seqClose(gds)
}
