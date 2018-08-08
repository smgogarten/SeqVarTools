library(Biobase)
library(GenomicRanges)

test_samplePairs1 <- function() {
  samples <- data.frame(sample.id=paste("samp", 1:10, sep=""),
                        subject.id=paste("subj", c(1,1,2,2,3,3,4,4,4,5), sep=""),
                        stringsAsFactors=FALSE)
  subj.df <- data.frame(sample1=paste("samp", c(1,3,5,7), sep=""),
                        sample2=paste("samp", c(2,4,6,8), sep=""),
                        row.names=paste("subj", 1:4, sep=""),
                        stringsAsFactors=FALSE)
  checkIdentical(subj.df, SeqVarTools:::.samplePairs1(samples))
}

test_samplePairs <- function() {
  samples <- data.frame(sample.id=paste("samp", 1:10, sep=""),
                        subject.id=paste("subj", c(1,1,2,2,3,3,4,4,4,5), sep=""),
                        stringsAsFactors=FALSE)
  subj.df <- data.frame(subject.id=paste("subj", 1:4, sep=""),
                        sample.id.1=paste("samp", c(1,3,5,7), sep=""),
                        sample.id.2=paste("samp", c(2,4,6,8), sep=""),
                        stringsAsFactors=FALSE)
  sp <- SeqVarTools:::.samplePairs(samples, all.pairs=FALSE)
  checkIdentical(subj.df, sp)
  sp1 <- suppressWarnings(SeqVarTools:::.samplePairs1(samples))
  checkIdentical(row.names(sp1), sp$subject.id)
  checkIdentical(sp1$sample1, sp$sample.id.1)
  checkIdentical(sp1$sample2, sp$sample.id.2)

  subj.df <- data.frame(subject.id=paste("subj", c(1:4,4,4), sep=""),
                        sample.id.1=paste("samp", c(1,3,5,7,7,8), sep=""),
                        sample.id.2=paste("samp", c(2,4,6,8,9,9), sep=""),
                        stringsAsFactors=FALSE)
  sp <- SeqVarTools:::.samplePairs(samples, all.pairs=TRUE)
  checkIdentical(subj.df, sp)
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

## test_duplicateDiscordance_apply <- function() {
##   gds <- SeqVarTools:::.testData()
##   sample.id <- seqGetData(gds, "sample.id")
##   samples <- data.frame(subject.id=c(rep(c("subj1", "subj2"), each=2), sample.id[5:length(sample.id)]),
##                         sample.id=sample.id,
##                         stringsAsFactors=FALSE)
##   var.id <- 101:110
  
##   seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(samples))
  
##   seqSetFilter(gds, variant.id=var.id, sample.id=sample.id[1:4])
  
##   disc <- duplicateDiscordance(seqData)
##   seqSetFilter(gds)
##   checkIdentical(disc,
##                  applyMethod(seqData, duplicateDiscordance,
##                              variant=var.id))
##   seqClose(gds)
## }



test_duplicateDiscordance <- function() {
  gds <- SeqVarTools:::.testData()
  sample.id <- seqGetData(gds, "sample.id")
  samples <- data.frame(subject.id=c(rep(c("subj1", "subj2"), each=2), sample.id[5:length(sample.id)]),
                        sample.id=sample.id,
                        stringsAsFactors=FALSE)
  var.id <- 101:110
  
  seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(samples))
  
  seqSetFilter(gds, variant.id=var.id, sample.id=sample.id[1:4], verbose=FALSE)
  disc1 <- SeqVarTools:::duplicateDiscordance1(seqData, verbose=FALSE)
  disc.subj <- duplicateDiscordance(seqData, by.variant=FALSE, all.pairs=FALSE, verbose=FALSE)
  disc.var <- duplicateDiscordance(seqData, by.variant=TRUE, all.pairs=FALSE, verbose=FALSE)

  checkEquals(disc.subj$subject.id, row.names(disc1$by.subject))
  checkEquals(disc.subj$sample.id.1, disc1$by.subject$sample1)
  checkEquals(disc.subj$sample.id.2, disc1$by.subject$sample2)
  checkEquals(disc.subj$n.variants, disc1$by.subject$num.var)
  checkEquals(disc.subj$n.concordant, disc1$by.subject$num.var - disc1$by.subject$num.discord)

  checkEquals(disc.var$variant.id, as.integer(rownames(disc1$by.variant)))
  checkEquals(disc.var$n.samples, disc1$by.variant$num.pair)
  checkEquals(disc.var$n.concordant, disc1$by.variant$num.pair - disc1$by.variant$num.discord)
  
  seqClose(gds)
}


test_duplicateDiscordance_iterator <- function() {
  gds <- SeqVarTools:::.testData()
  sample.id <- seqGetData(gds, "sample.id")
  samples <- data.frame(subject.id=c(rep(c("subj1", "subj2"), each=3), sample.id[7:length(sample.id)]),
                        sample.id=sample.id,
                        stringsAsFactors=FALSE)
  
  seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(samples))
  disc.subj <- duplicateDiscordance(seqData, by.variant=FALSE, all.pairs=TRUE, verbose=FALSE)
  disc.var <- duplicateDiscordance(seqData, by.variant=TRUE, all.pairs=TRUE, verbose=FALSE)

  it <- SeqVarBlockIterator(seqData, variantBlock=500, verbose=FALSE)
  disc.subj.it <- duplicateDiscordance(it, by.variant=FALSE, all.pairs=TRUE, verbose=FALSE)
  resetIterator(it, verbose=FALSE)
  disc.var.it <- duplicateDiscordance(it, by.variant=TRUE, all.pairs=TRUE, verbose=FALSE)

  checkEquals(disc.subj, disc.subj.it)
  checkEquals(disc.var, disc.var.it)

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

test_getMatchesHetHom <- function (){
  df <- .getTestGeno()
  truth <- c(T, F, T, F, F, T, F, F, T, F, T, F, F, F, F, F)
  checkEquals(SeqVarTools:::.getMatchesHetHom(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getMatchesHetHom(df$geno1, df$geno2), SeqVarTools:::.getMatchesHetHom(df$geno2, df$geno1))
}

test_getNonMissing <- function (){
  df <- .getTestGeno()
  truth <- c(T, T, T, F, T, T, T, F, T, T, T, F, F, F, F, F)
  checkEquals(SeqVarTools:::.getNonMissing(df$geno1, df$geno2), truth)
  checkIdentical(SeqVarTools:::.getNonMissing(df$geno1, df$geno2), SeqVarTools:::.getNonMissing(df$geno2, df$geno1))
}

test_getGentoypeClass <- function() {
  vals <- c(NA, 0, 1, 2)
  class <- SeqVarTools:::.getGenotypeClass(vals)
  checkEquals(class, c("miss", "alt", "het", "ref"))

}

test_matchVariants <- function() {
  
  gr1 <- GRanges(seqnames=c("1", "1", "1", "1", "2", "2", "4"),
                 ranges=IRanges(start=c(1, 1, 2, 3, 1, 3, 1),
                                end  =c(1, 2, 2, 3, 1, 3, 1))
  )
  gr1$ref <- rep("A", length(gr1))
  gr1$alt <- rep("B", length(gr1))
  gr1$alt[2] <- "BC"
  gr1$snv <- rep(TRUE, length(gr1))
  gr1$snv[2] <- FALSE
  gr1$variant.id <- 1:length(gr1) + 2
  
  gr2 <- GRanges(seqnames=c("1", "1", "1", "1", "1", "2", "2", "3", "4"),
                 ranges=IRanges(start=c(1, 1, 1, 2, 2, 1, 3, 1, 1),
                                end = c(1, 1, 2, 2, 2, 1, 3, 1, 1)))
  gr2$ref <- rep("A", length(gr2))
  gr2$alt <- rep("B", length(gr2))
  gr2$alt[3] <- "BC"
  gr2$alt[4] <- "C"
  gr2$ref[9] <- "B"
  gr2$alt[9] <- "A"
  gr2$snv <- rep(TRUE, length(gr2))
  gr2$snv[3] <- FALSE
  gr2$variant.id <- 1:length(gr2)
  
  # match on position only, duplicates alowed
  overlaps <- SeqVarTools:::.matchVariants(gr1, gr2, match.on="position")  
  checkEquals(overlaps$variant.id.1, c(3, 3, 5, 5, 7, 8, 9))
  checkEquals(overlaps$variant.id.2, c(1, 2, 4, 5, 6, 7, 9))

  # match on position, no duplicates
  overlapsNoDups <- SeqVarTools:::.matchVariants(gr1, gr2, match.on="position", allowOverlaps=FALSE)
  checkEquals(overlapsNoDups$variant.id.1, c(3, 5, 7, 8, 9))
  checkEquals(overlapsNoDups$variant.id.2, c(1, 4, 6, 7, 9))
  
  # match on position and alleles  
  overlaps <- SeqVarTools:::.matchVariants(gr1, gr2, match.on="alleles")
  checkEquals(overlaps$variant.id.1, c(3, 3, 5, 7, 8, 9))
  checkEquals(overlaps$variant.id.2, c(1, 2, 5, 6, 7, 9))
  checkEquals(overlaps$recode, c(rep(FALSE, 5), TRUE))

  overlapsNoDups <- SeqVarTools:::.matchVariants(gr1, gr2, match.on="alleles", allowOverlaps=FALSE)
  checkEquals(overlapsNoDups$variant.id.1, c(3, 5, 7, 8, 9))
  checkEquals(overlapsNoDups$variant.id.2, c(1, 5, 6, 7, 9))
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
  res <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles")
  res.var <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", by.variant=TRUE)
  
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
  
  res <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles") 
  res.var <- duplicateDiscordance(seqData1, seqData2, by.variant=TRUE, match.variants.on="alleles") 
  
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



test_duplicateDiscordanceAcrossDatasets_hethom <- function() {
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
  res <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", discordance.type = "hethom")
  tmp <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", discordance.type = "genotype")
  checkEquals(res$n.variants, tmp$n.variants)
  checkEquals(res$n.concordant, tmp$n.concordant)
  
  res.var <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", by.variant=TRUE, discordance.type="hethom")
  tmp <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", by.variant=TRUE)
  checkEquals(res.var$n.samples, tmp$n.samples)
  checkEquals(res.var$n.concordant, tmp$n.concordant)
  
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
  
  seqSetFilter(seqData1)
  seqSetFilter(seqData2)
  
  # change subjectID for seqData2
  #samples2 <- sampleData(seqData2)
  samples2$subject.id[2:1] <- samples2$subject.id[1:2]
  seqData2 <- SeqVarData(gds2, sampleData=AnnotatedDataFrame(samples2))
  
  
  seqSetFilter(seqData1, sample.id=samples1$sample.id[1:2], variant.id=seqGetData(seqData1, "variant.id")[1:50])
  seqSetFilter(seqData2, sample.id=samples2$sample.id[1:2], variant.id=seqGetData(seqData2, "variant.id")[25:75])
  
  res <- duplicateDiscordance(seqData1, seqData2, match.variants.on="alleles", discordance.type="hethom") 
  res.var <- duplicateDiscordance(seqData1, seqData2, by.variant=TRUE, match.variants.on="alleles", discordance.type="hethom") 
  
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
  dos <- abs(dos - 1)
  match <- dos[1 ,] == dos[2, ]
  checkEquals(res$n.concordant[1], sum(match))
  
  seqClose(seqData1)
  seqClose(seqData2)
  
  unlink(tmpfile)
  
}
