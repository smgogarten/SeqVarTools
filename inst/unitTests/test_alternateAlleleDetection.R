library(Biobase)

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
  gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
  
  filename <- seqExampleFileName("gds")
  tmpfile <- tempfile()
  file.copy(filename, tmpfile)
  
  gds1 <- seqOpen(filename)
  gds2 <- seqOpen(tmpfile)
  
  # add sample data
  samples1 <- data.frame(sample.id=seqGetData(gds1, "sample.id"), stringsAsFactors=F)
  samples1$subject.id <- samples1$sample.id
  seqData1 <- SeqVarData(gds1, sampleData=AnnotatedDataFrame(samples1))
  
  samples2 <- data.frame(sample.id=seqGetData(gds2, "sample.id"), stringsAsFactors=F)
  samples2$subject.id <- samples2$sample.id
  seqData2 <- SeqVarData(gds2, sampleData=AnnotatedDataFrame(samples2))
  
  seqSetFilter(seqData1)
  seqSetFilter(seqData2)
  
  seqSetFilter(seqData1, sample.id=seqGetData(seqData1, "sample.id")[1:3])
  seqSetFilter(seqData2, sample.id=seqGetData(seqData2, "sample.id")[2:4])
  
  # check filters
  filt.1 <- seqGetFilter(seqData1)
  filt.2 <- seqGetFilter(seqData2)
  
  # makes sure it runs
  res <- alternateAlleleDetection(seqData1, seqData2) 
  
  # check filters
  checkEquals(seqGetFilter(seqData1), filt.1)
  checkEquals(seqGetFilter(seqData2), filt.2)
  
  # check data
  checkTrue(all(res$false.pos == 0))
  checkTrue(all(res$false.neg == 0))
  
  samps <- intersect(seqGetData(seqData1, "sample.id"), seqGetData(seqData2, "sample.id"))
  
  seqSetFilter(seqData1, sample.id=samps, variant.id=res$variant.id.1)
  refdos <- refDosage(seqData1)
  # number of samples with nonmissing data
  checkEquals(res$n.samples, colSums(!is.na(refdos)), checkNames=FALSE)
  checkEquals(res$true.neg, colSums(refdos, na.rm=T), checkNames=FALSE)
  checkEquals(res$true.pos, colSums(2-refdos, na.rm=T), checkNames=FALSE)

  # change the mapping
  samples2$subject.id[2:1] <- samples2$subject.id[1:2]
  seqData2 <- SeqVarData(gds2, sampleData=AnnotatedDataFrame(samples2))
  
  
  seqSetFilter(seqData1)
  seqSetFilter(seqData2)
  
  seqSetFilter(seqData1, sample.id=samples1$sample.id[1:2], variant.id=seqGetData(seqData1, "variant.id")[1:50])
  seqSetFilter(seqData2, sample.id=samples2$sample.id[1:2], variant.id=seqGetData(seqData2, "variant.id")[25:75])
  
  res <- alternateAlleleDetection(seqData1, seqData2) 
  
  checkEquals(res$variant.id.1, intersect(seqGetData(seqData1, "variant.id"), seqGetData(seqData2, "variant.id")))
  
  # set filter to only read overlaps
  seqSetFilter(seqData1, variant.id=res$variant.id.1)
  dos <- refDosage(seqData1)
  
  class1 <- SeqVarTools:::.getGenotypeClass(dos[1, ])
  class2 <- SeqVarTools:::.getGenotypeClass(dos[2, ])
  
  checkEquals(res$true.pos, 2*SeqVarTools:::.truePos(class1, class2))
  checkEquals(res$true.neg, 2*SeqVarTools:::.trueNeg(class1, class2))
  checkEquals(res$false.pos, SeqVarTools:::.falsePos(class1, class2) + SeqVarTools:::.falsePos(class2, class1))
  checkEquals(res$false.neg, SeqVarTools:::.falseNeg(class1, class2) + SeqVarTools:::.falseNeg(class2, class1))

  seqClose(seqData1)
  seqClose(seqData2)
  
  unlink(tmpfile)
  
}
