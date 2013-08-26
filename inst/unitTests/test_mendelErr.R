test_triosFromPedigree <- function() {
  ped <- data.frame(family=c(rep("a",4), rep("b",5), rep("c",2)),
                    individ=c("a1","a2","am","af", "b1","bm","bf","bfm","bff", "c1","cm"),
                    father=c("af","af",0,0, "bf",0,"bff",0,0, 0,0),
                    mother=c("am","am",0,0, "bm",0,"bfm",0,0, "cm",0),
                    sample.id=c(11:14, 21:25, 31:32),
                    stringsAsFactors=FALSE)
  trios <- list(c(child=11, mother=13, father=14),
                c(child=12, mother=13, father=14),
                c(child=21, mother=22, father=23),
                c(child=23, mother=24, father=25),
                c(child=31, mother=32))
  checkEquals(trios, SeqVarTools:::.triosFromPedigree(ped))
}


test_alleleMatch <- function() {
  geno <- array(dim=c(2,2,5),
                dimnames=list(NULL,c("child","parent"),NULL))
  geno[1,1,] <- c(0,0,0,0,0)
  geno[2,1,] <- c(1,1,1,1,1)
  geno[1,2,] <- c(0,0,0,1,1)
  geno[2,2,] <- c(0,0,0,0,1)
  match1 <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
  match2 <- c(FALSE,FALSE,FALSE,TRUE,TRUE)
  checkIdentical(match1, SeqVarTools:::.alleleMatch(geno,1,"parent"))
  checkIdentical(match2, SeqVarTools:::.alleleMatch(geno,2,"parent"))
}


test_autosomeErr <- function() {
  geno <- array(dim=c(2,3,9),
                dimnames=list(NULL, sample=c("child","mother","father"), NULL))
  ## father  mother  child
  ## AB      AB      AB     1 - no err
  ## AB      AB      AA     2 - no err
  ## AB      AB      BB     3 - no err
  ## AA      BB      AB     4 - no err
  ## AA      BB      AA     5 - err
  ## AA      BB      BB     6 - err
  ## BB      AA      AB     7 - no err
  ## BB      AA      AA     8 - err
  ## BB      AA      BB     9 - err
  geno[1,"father",] <- c(0,0,0,0,0,0,1,1,1)
  geno[2,"father",] <- c(1,1,1,0,0,0,1,1,1)
  geno[1,"mother",] <- c(0,0,0,1,1,1,0,0,0)
  geno[2,"mother",] <- c(1,1,1,1,1,1,0,0,0)
  geno[1,"child",]  <- c(0,0,1,0,0,1,0,0,1)
  geno[2,"child",]  <- c(1,0,1,1,0,1,1,0,1) 
  err <- c(rep(FALSE,3),rep(c(FALSE,TRUE,TRUE),2))
  checkIdentical(err, SeqVarTools:::.autosomeErr(geno))
  err.mother <- c(rep(FALSE,3),FALSE,TRUE,FALSE,FALSE,FALSE,TRUE)
  checkIdentical(err.mother, SeqVarTools:::.autosomeErr(geno[,c("child","mother"),]))
  err.father <- c(rep(FALSE,3),FALSE,FALSE,TRUE,FALSE,TRUE,FALSE)
  checkIdentical(err.father, SeqVarTools:::.autosomeErr(geno[,c("child","father"),]))
                                          
  geno <- array(dim=c(2,3,6),
                dimnames=list(NULL, sample=c("child","mother","father"), NULL))
  ## father  mother  child
  ## AA      AB      AB     1 - no err
  ## AA      AB      AA     2 - no err
  ## AA      AB      BB     3 - err
  ## AB      AA      AB     4 - no err
  ## AB      AA      AA     5 - no err
  ## AB      AA      BB     6 - err
  geno[1,"father",] <- c(0,0,0,0,0,0)
  geno[2,"father",] <- c(0,0,0,1,1,1)
  geno[1,"mother",] <- c(0,0,0,0,0,0)
  geno[2,"mother",] <- c(1,1,1,0,0,0)
  geno[1,"child",]  <- c(0,0,1,0,0,1)
  geno[2,"child",]  <- c(1,0,1,1,0,1) 
  err <- rep(c(FALSE,FALSE,TRUE),2)
  checkIdentical(err, SeqVarTools:::.autosomeErr(geno))
  err.mother <- c(rep(FALSE,3),FALSE,FALSE,TRUE)
  checkIdentical(err.mother, SeqVarTools:::.autosomeErr(geno[,c("child","mother"),]))
  err.father <- c(FALSE,FALSE,TRUE,rep(FALSE,3))
  checkIdentical(err.father, SeqVarTools:::.autosomeErr(geno[,c("child","father"),]))
  
  geno <- array(dim=c(2,3,6),
                dimnames=list(NULL, sample=c("child","mother","father"), NULL))
  ## father  mother  child
  ## AA      **      AB     1 - no err
  ## AA      **      AA     2 - no err
  ## AA      **      BB     3 - err
  ## **      AA      AB     4 - no err
  ## **      AA      AA     5 - no err
  ## **      AA      BB     6 - err
  geno[1,"father",] <- c(0,0,0,NA,NA,NA)
  geno[2,"father",] <- c(0,0,0,NA,NA,NA)
  geno[1,"mother",] <- c(NA,NA,NA,0,0,0)
  geno[2,"mother",] <- c(NA,NA,NA,0,0,0)
  geno[1,"child",]  <- c(0,0,1,0,0,1)
  geno[2,"child",]  <- c(1,0,1,1,0,1) 
  err <- rep(c(FALSE,FALSE,TRUE),2)
  checkIdentical(err, SeqVarTools:::.autosomeErr(geno))
}

test_xMaleErr <- function() {
  geno <- array(dim=c(2,2,6),
                dimnames=list(NULL, sample=c("child","mother"), NULL))
  ## father  mother  child
  ## **      AA      AB    1 - no err
  ## **      AA      AA    2 - no err
  ## **      AA      BB    3 - err
  ## **      BB      AB    4 - no err
  ## **      BB      AA    5 - err
  ## **      BB      BB    6 - no err
  geno[1,"mother",] <- c(0,0,0,1,1,1)
  geno[2,"mother",] <- c(0,0,0,1,1,1)
  geno[1,"child",]  <- c(0,0,1,0,0,1)
  geno[2,"child",]  <- c(1,0,1,1,0,1)
  err <- c(FALSE,FALSE,TRUE,FALSE,TRUE,FALSE)
  checkIdentical(err, SeqVarTools:::.xMaleErr(geno))
}

test_yMaleErr <- function() {
  geno <- array(dim=c(2,2,6),
                dimnames=list(NULL, sample=c("child","father"), NULL))
  ## father  mother  child
  ## AA      **      AB    1 - no err
  ## AA      **      AA    2 - no err
  ## AA      **      BB    3 - err
  ## BB      **      AB    4 - no err
  ## BB      **      AA    5 - err
  ## BB      **      BB    6 - no err
  geno[1,"father",] <- c(0,0,0,1,1,1)
  geno[2,"father",] <- c(0,0,0,1,1,1)
  geno[1,"child",]  <- c(0,0,1,0,0,1)
  geno[2,"child",]  <- c(1,0,1,1,0,1)
  err <- c(FALSE,FALSE,TRUE,FALSE,TRUE,FALSE)
  checkIdentical(err, SeqVarTools:::.yMaleErr(geno))
}

test_mendelErr <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  data(pedigree)
  checkEquals(list(c(child="NA12878", mother="NA12892", father="NA12891")),
              SeqVarTools:::.triosFromPedigree(pedigree))
  err <- mendelErr(gds, pedigree)
  checkEquals(rep(0, length(seqGetData(gds, "variant.id"))),
              err$by.variant)
  checkEquals(c("NA12878"=0), err$by.trio)
  seqClose(gds)
}

test_mendelErr_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  data(pedigree)
  var.id <- 101:110
  seqSetFilter(gds, variant.id=var.id)
  err <- mendelErr(gds, pedigree)
  seqSetFilter(gds)
  checkIdentical(err,
                 applyMethod(gds, mendelErr, variant=var.id, pedigree=pedigree))
  seqClose(gds)
}
