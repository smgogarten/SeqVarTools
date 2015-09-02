test_hwe <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  pv <- hwe(gds, use.names=TRUE)
  checkIdentical(as.character(seqGetData(gds, "variant.id")),
                 names(pv))
  
  mono <- alleleFrequency(gds) %in% c(0,1)
  checkEquals(mono | nAlleles(gds) != 2, is.na(pv), checkNames=FALSE)
  seqClose(gds)
}

test_hwe_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  pv <- hwe(gds)
  seqSetFilter(gds)
  checkIdentical(pv,
                 applyMethod(gds, hwe, variant=var.id, sample=samp.id))
  seqClose(gds)
}

