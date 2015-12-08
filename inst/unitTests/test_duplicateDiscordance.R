library(Biobase)
library(GenomicRanges)

test_samplePairs <- function() {
  samples <- data.frame(sample.id=paste("samp", 1:10, sep=""),
                        subject.id=paste("subj", c(1,1,2,2,3,3,4,4,4,5), sep=""),
                        stringsAsFactors=FALSE)
  subj.df <- data.frame(sample1=paste("samp", c(1,3,5,7), sep=""),
                        sample2=paste("samp", c(2,4,6,8), sep=""),
                        row.names=paste("subj", 1:4, sep=""),
                        stringsAsFactors=FALSE)
  checkIdentical(subj.df, SeqVarTools:::.samplePairs(samples))
}

test_genoMatch <- function() {
  geno <- array(dim=c(2,2,5))
  geno[1,1,] <- c(0,0,0,0,0)
  geno[2,1,] <- c(1,1,1,1,1)
  geno[1,2,] <- c(0,1,0,1,NA)
  geno[2,2,] <- c(1,0,0,1,NA)
  match.phase <- c(TRUE,FALSE,FALSE,FALSE,NA)
  match.unphase <- c(TRUE,TRUE,FALSE,FALSE,NA)
  checkIdentical(match.phase, SeqVarTools:::.genoMatch(geno, check.phase=TRUE))
  checkIdentical(match.unphase, SeqVarTools:::.genoMatch(geno, check.phase=FALSE))
}

test_duplicateDiscordance_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  sample.id <- seqGetData(gds, "sample.id")
  samples <- data.frame(subject.id=c(rep(c("subj1", "subj2"), each=2), sample.id[5:length(sample.id)]),
                        sample.id=sample.id,
                        stringsAsFactors=FALSE)
  var.id <- 101:110
  
  seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(samples))
  
  seqSetFilter(gds, variant.id=var.id, sample.id=sample.id[1:4])
  
  disc <- duplicateDiscordance(seqData)
  seqSetFilter(gds)
  checkIdentical(disc,
                 applyMethod(seqData, duplicateDiscordance,
                             variant=var.id))
  seqClose(gds)
}



## tests for duplicate discordance on two gds files

.getTestGeno <- function() {
  vals <- c("ref", "het", "alt", "miss")
  df <- expand.grid(geno1=vals, geno2=vals, stringsAsFactors=F)
  
  #n.random <- 20
  #df2 <- data.frame(geno1=sample(vals, n.random, replace=T), geno2=sample(vals, n.random, replace=T), stringsAsFactors=F)
  #df <- rbind(df, df2)
  df  
}

test_getMatchesConc <- function(){
  df <- .getTestGeno()
  truth <- c(T, F, F, F, F, T, F, F, F, F, T, F, F, F, F, F)
  checkEquals(SeqVarTools:::.getMatchesConc(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesConc(df$geno1, df$geno2), SeqVarTools:::.getMatchesConc(df$geno2, df$geno1))
}


test_getAlt <- function(){
  df <- .getTestGeno()
  truth <- c(F, T, T, F, T, T, T, F, T, T, T, F, F, F, F, F)
  
  checkEquals(SeqVarTools:::.getAlt(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getAlt(df$geno1, df$geno2), SeqVarTools:::.getAlt(df$geno2, df$geno1))  
}


test_getMatchesAltConc <- function(){
  df <- .getTestGeno()
  truth <- c(F, F, F, F, F, T, F, F, F, F, T, F, F, F, F, F)
  
  checkEquals(SeqVarTools:::.getMatchesAltConc(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesAltConc(df$geno1, df$geno2), SeqVarTools:::.getMatchesAltConc(df$geno2, df$geno1))
}


test_getMatchesHetRef <- function(){
  df <- .getTestGeno()
  
  truth <- c(F, T, F, F, T, F, F, F, F, F, F, F, F, F, F, F)

  checkEquals(SeqVarTools:::.getMatchesHetRef(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesHetRef(df$geno1, df$geno2), SeqVarTools:::.getMatchesHetRef(df$geno2, df$geno1))
}

test_getMatchesHetAlt <- function(){
  df <- .getTestGeno()
  truth <- c(F, F, F, F, F, F, T, F, F, T, F, F, F, F, F, F)
    
  checkEquals(SeqVarTools:::.getMatchesHetAlt(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesHetAlt(df$geno1, df$geno2), SeqVarTools:::.getMatchesHetAlt(df$geno2, df$geno1))
}


test_getMatchesRefAlt <- function() {
  df <- .getTestGeno()
  truth <- c(F, F, T, F, F, F, F, F, T, F, F, F, F, F, F, F)
  checkEquals(SeqVarTools:::.getMatchesRefAlt(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesRefAlt(df$geno1, df$geno2), SeqVarTools:::.getMatchesRefAlt(df$geno2, df$geno1))
}

test_getGentoypeClass <- function() {
  
  vals <- c(NA, 0, 1, 2)
  class <- SeqVarTools:::.getGenotypeClass(vals)
  checkEquals(class, c("miss", "alt", "het", "ref"))

}

test_matchVariants <- function() {
  
  gr1 <- GRanges(seqnames=c("1", "1", "1", "1", "2", "2"),
                 ranges=IRanges(start=c(1, 1, 2, 3, 1, 3),
                                end  =c(1, 2, 2, 3, 1, 3))
  )
  gr1$ref <- rep("A", length(gr1))
  gr1$alt <- rep("B", length(gr1))
  gr1$alt[2] <- "BC"
  gr1$snv <- rep(TRUE, length(gr1))
  gr1$snv[2] <- FALSE
  gr1$variant.id <- 1:length(gr1) + 2
  
  gr2 <- GRanges(seqnames=c("1", "1", "1", "1", "1", "2", "2", "3"),
                 ranges=IRanges(start=c(1, 1, 1, 2, 2, 1, 3, 1),
                                end = c(1, 1, 2, 2, 2, 1, 3, 1)))
  gr2$ref <- rep("A", length(gr2))
  gr2$alt <- rep("B", length(gr2))
  gr2$alt[3] <- "BC"
  gr2$alt[4] <- "C"
  gr2$snv <- rep(TRUE, length(gr2))
  gr2$snv[3] <- FALSE
  gr2$variant.id <- 1:length(gr2)
  
  overlaps <- SeqVarTools:::.matchVariants(gr1, gr2)  
  
  checkEquals(overlaps$variant.id.1, c(3, 3, 5, 7, 8))
  checkEquals(overlaps$variant.id.2, c(1, 2, 5, 6, 7))
  
  overlapsNoDups <- SeqVarTools:::.matchVariants(gr1, gr2, allowOverlaps=FALSE)
  checkEquals(overlapsNoDups$variant.id.1, c(3, 5, 7, 8))
  checkEquals(overlapsNoDups$variant.id.2, c(1, 5, 6, 7))
}

test_matchSamples <- function() {
  
  # no duplicates at this point
  samp1 <- data.frame(subject.id=c("A", "A", "B", "C", "D"), stringsAsFactors=F)
  samp1$sample.id <- letters[1:nrow(samp1)]                    
  
  samp2 <- data.frame(subject.id=c("A", "B", "D", "D", "E"), stringsAsFactors=F)
  samp2$sample.id <- letters[1:nrow(samp2)]
  
  chk <- SeqVarTools:::.matchSamples(samp1, samp2)
  chk$chk <- paste(chk$subject.id, chk$sample.id.1, chk$sample.id.2)
  checkTrue(setequal(chk$chk, c("A a a", "A b a", "B c b", "D e c", "D e d")))
  checkTrue(setequal(chk$subject.id, intersect(samp1$subject.id, samp2$subject.id)))
  checkEquals(nrow(chk), 5)
  
}


test_duplicateDiscordanceAcrossDatasets <- function() {
  
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
  res <- duplicateDiscordance(seqData1, seqData2)
  res.var <- duplicateDiscordance(seqData1, seqData2, by.variant=TRUE)
  
  # check filters
  checkEquals(seqGetFilter(seqData1), filt.1)
  checkEquals(seqGetFilter(seqData2), filt.2)
  
  # check data for samples
  checkEquals(seqGetData(seqData1, "sample.id")[2:3], res$sample.id.1)
  checkEquals(seqGetData(seqData2, "sample.id")[1:2], res$sample.id.2)
  # and for variants
  checkEquals(seqGetData(seqData1, "variant.id")[isSNV(seqData1)], res.var$variant.id.1)
  checkEquals(seqGetData(seqData2, "variant.id")[isSNV(seqData2)], res.var$variant.id.2)
  
  checkEquals(res$n.concordant, res$n.variants)
  checkEquals(res$n.alt, res$n.alt.conc)

  checkEquals(res.var$n.concordant, res.var$n.samples)
  checkEquals(res.var$n.alt, res.var$n.alt.conc)
  
  seqSetFilter(seqData1)
  seqSetFilter(seqData2)

    # change subjectID for seqData2
  #samples2 <- sampleData(seqData2)
  samples2$subject.id[2:1] <- samples2$subject.id[1:2]
  seqData2 <- SeqVarData(gds2, sampleData=AnnotatedDataFrame(samples2))
  
  
  seqSetFilter(seqData1, sample.id=samples1$sample.id[1:2], variant.id=seqGetData(seqData1, "variant.id")[1:50])
  seqSetFilter(seqData2, sample.id=samples2$sample.id[1:2], variant.id=seqGetData(seqData2, "variant.id")[25:75])
  
  res <- duplicateDiscordance(seqData1, seqData2) 
  res.var <- duplicateDiscordance(seqData1, seqData2, by.variant=TRUE) 
  
  # set filter to only read overlaps
  seqSetFilter(seqData1, variant.id=intersect(seqGetData(seqData1, "variant.id"), seqGetData(gds2, "variant.id")))
  
  # read dosages
  dos <- refDosage(seqData1)
  
  # get the granges object
  gr <- SeqVarTools:::.getGRanges(seqData1)
  
  # snvs only
  dos <- dos[, gr$snv]
  
  # non missing
  nonmiss <- !is.na(colSums(dos))
  checkEquals(res$n.variants[1], sum(nonmiss))
  
  dos <- dos[, nonmiss]
  match <- dos[1 ,] == dos[2, ]
  checkEquals(res$n.concordant[1], sum(match))
  
  # minor allele matches
  minor <- dos[1, ] < 2 | dos[2, ] < 2
  
  checkEquals(res$n.alt.conc[1], sum(minor & match))
  
  seqClose(seqData1)
  seqClose(seqData2)
  
  unlink(tmpfile)
  
}
