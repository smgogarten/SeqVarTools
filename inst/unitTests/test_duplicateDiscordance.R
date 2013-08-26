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
  samples <- data.frame(subject.id=rep(c("subj1", "subj2"), each=2),
                        sample.id=sample.id[1:4],
                        stringsAsFactors=FALSE)
  var.id <- 101:110
  seqSetFilter(gds, variant.id=var.id)
  disc <- duplicateDiscordance(gds, samples)
  seqSetFilter(gds)
  checkIdentical(disc,
                 applyMethod(gds, duplicateDiscordance,
                             variant=var.id, samples=samples))
  seqClose(gds)
}
