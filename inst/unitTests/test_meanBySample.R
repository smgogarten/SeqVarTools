test_meanBySample <- function() {
  gds <- SeqVarTools:::.testData()
  dp <- getVariableLengthData(gds, "annotation/format/DP", use.names=TRUE)
  checkEquals(rowMeans(dp, na.rm=TRUE),
              meanBySample(gds, "annotation/format/DP", use.names=TRUE))
  seqClose(gds)
}

test_meanBySample_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  mn <- meanBySample(gds, "annotation/format/DP")
  seqSetFilter(gds)
  checkIdentical(mn,
                 applyMethod(gds, meanBySample, variant=var.id, sample=samp.id,
                             var.name="annotation/format/DP"))
  seqClose(gds)
}
