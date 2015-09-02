
test_posNeg <- function() {
  TP <- matrix(c(2,1,0,0,
                 1,1,0,0,
                 0,0,0,0,
                 0,0,0,0), nrow=4, byrow=TRUE)
  TN <- matrix(c(0,0,0,0,
                 0,1,1,0,
                 0,1,2,0,
                 0,0,0,0), nrow=4, byrow=TRUE)
  FP <- matrix(c(0,0,0,0,
                 1,0,0,0,
                 2,1,0,0,
                 0,0,0,0), nrow=4, byrow=TRUE)
  FN <- matrix(c(0,1,2,0,
                 0,0,1,0,
                 0,0,0,0,
                 0,0,0,0), nrow=4, byrow=TRUE)
  
  geno1 <- matrix(c("alt", "het", "ref", "miss"),
                  nrow=4, ncol=4, byrow=TRUE)
  geno2 <-  matrix(c("alt", "het", "ref", "miss"),
                   nrow=4, ncol=4)
  
  checkIdentical(TP, SeqVarTools:::.truePos(geno1, geno2))
  checkIdentical(TN, SeqVarTools:::.trueNeg(geno1, geno2))
  checkIdentical(FP, SeqVarTools:::.falsePos(geno1, geno2))
  checkIdentical(FN, SeqVarTools:::.falseNeg(geno1, geno2))
  
}



test_alternateAlleleDetection <- function() {
  
  filename <- seqExampleFileName("gds")
  tmpfile <- tempfile()
  file.copy(filename, tmpfile)
  
  gds1 <- seqOpen(filename)
  gds2 <- seqOpen(tmpfile)
  
  seqSetFilter(gds1)
  seqSetFilter(gds2)
  
  seqSetFilter(gds1, sample.id=seqGetData(gds1, "sample.id")[1:3])
  seqSetFilter(gds2, sample.id=seqGetData(gds2, "sample.id")[2:4])
  
  # check filters
  filt.1 <- seqGetFilter(gds1)
  filt.2 <- seqGetFilter(gds2)
  
  # makes sure it runs
  res <- alternateAlleleDetection(gds1, gds2) 
  
  # check filters
  checkEquals(seqGetFilter(gds1), filt.1)
  checkEquals(seqGetFilter(gds2), filt.2)
  
  # check data
  checkTrue(all(res$false.pos == 0))
  checkTrue(all(res$false.neg == 0))
  
  samps <- intersect(seqGetData(gds1, "sample.id"), seqGetData(gds2, "sample.id"))
  
  seqSetFilter(gds1, sample.id=samps, variant.id=res$variant.id.1)
  refdos <- refDosage(gds1)
  # number of samples with nonmissing data
  checkEquals(res$n.samples, colSums(!is.na(refdos)), checkNames=FALSE)
  checkEquals(res$true.neg, colSums(refdos, na.rm=T), checkNames=FALSE)
  checkEquals(res$true.pos, colSums(2-refdos, na.rm=T), checkNames=FALSE)
  closefn.gds(gds2)
  
  # change a sample id
  tmp <- openfn.gds(tmpfile, readonly=FALSE)  
  samps <- read.gdsn(index.gdsn(tmp, "sample.id"))
  samps[1:2] <- samps[2:1]
  delete.gdsn(index.gdsn(tmp, "sample.id"))
  add.gdsn(tmp, "sample.id", val=samps, compress="ZIP")
  closefn.gds(tmp)
  
  gds2 <- seqOpen(tmpfile)  
  
  seqSetFilter(gds1)
  seqSetFilter(gds2)
  
  seqSetFilter(gds1, sample.id=samps[1:2], variant.id=seqGetData(gds1, "variant.id")[1:50])
  seqSetFilter(gds2, sample.id=samps[1:2], variant.id=seqGetData(gds2, "variant.id")[25:75])
  
  res <- alternateAlleleDetection(gds1, gds2) 
  
  checkEquals(res$variant.id.1, intersect(seqGetData(gds1, "variant.id"), seqGetData(gds2, "variant.id")))
  
  # set filter to only read overlaps
  seqSetFilter(gds1, variant.id=res$variant.id.1)
  dos <- refDosage(gds1)
  
  class1 <- SeqVarTools:::.getGenotypeClass(dos[1, ])
  class2 <- SeqVarTools:::.getGenotypeClass(dos[2, ])
  
  checkEquals(res$true.pos, 2*SeqVarTools:::.truePos(class1, class2))
  checkEquals(res$true.neg, 2*SeqVarTools:::.trueNeg(class1, class2))
  checkEquals(res$false.pos, SeqVarTools:::.falsePos(class1, class2) + SeqVarTools:::.falsePos(class2, class1))
  checkEquals(res$false.neg, SeqVarTools:::.falseNeg(class1, class2) + SeqVarTools:::.falseNeg(class2, class1))

  closefn.gds(gds1)
  closefn.gds(gds2)
  
  unlink(tmpfile)
  
}
