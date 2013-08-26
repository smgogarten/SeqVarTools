test_applyNames <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  samp.id <- as.character(seqGetData(gds, "sample.id"))
  var.id <- as.character(seqGetData(gds, "variant.id"))
  
  var <- matrix(nrow=length(samp.id), ncol=length(var.id),
                dimnames=list(sample=NULL, variant=NULL))
  checkIdentical(list(sample=samp.id, variant=var.id),
                 dimnames(SeqVarTools:::.applyNames(gds, var)))
  
  var <- matrix(nrow=length(var.id), ncol=length(samp.id),
                dimnames=list(variant=NULL, sample=NULL))
  checkIdentical(list(variant=var.id, sample=samp.id),
                 dimnames(SeqVarTools:::.applyNames(gds, var)))

  var <- array(dim=c(2, length(samp.id), length(var.id)),
               dimnames=list(NULL, sample=NULL, variant=NULL))
  checkIdentical(list(NULL, sample=samp.id, variant=var.id),
                 dimnames(SeqVarTools:::.applyNames(gds, var)))
  
  var <- matrix(nrow=10, ncol=10)
  checkIdentical(NULL, dimnames(SeqVarTools:::.applyNames(gds, var)))
  seqClose(gds)
}


test_isVariant <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  var <- (!is.na(geno[1,,]) & geno[1,,] > 0) | (!is.na(geno[2,,]) & geno[2,,] > 0)
  dimnames(var) <- list(sample=NULL, variant=NULL)
  checkIdentical(var, isVariant(gds))
  seqClose(gds)
}

test_isVariant_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  var <- isVariant(gds)
  seqSetFilter(gds)
  checkIdentical(var,
                 applyMethod(gds, isVariant, variant=var.id, sample=samp.id))
  seqClose(gds)
}


test_getGenotype <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  gc <- matrix(paste(geno[1,,], geno[2,,], sep="/"),
               nrow=dim(geno)[2], ncol=dim(geno)[3],
               dimnames=list(sample=NULL, variant=NULL))
  gc[gc == "NA/NA"] <- NA
  checkIdentical(gc, getGenotype(gds, use.names=FALSE))
  seqClose(gds)
}

test_getGenotype_phased <- function() {
  gds <- seqOpen(system.file("extdata", "gl_chr1.gds", package="SeqVarTools"))
  geno <- seqGetData(gds, "genotype")
  gc <- matrix(paste(geno[1,,], geno[2,,], sep="|"),
               nrow=dim(geno)[2], ncol=dim(geno)[3],
               dimnames=list(sample=NULL, variant=NULL))
  gc[gc == "NA/NA"] <- NA
  checkIdentical(gc, getGenotype(gds, use.names=FALSE))
  seqClose(gds)
}

test_getGenotype_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  geno <- getGenotype(gds)
  seqSetFilter(gds)
  checkIdentical(geno,
                 applyMethod(gds, getGenotype, variant=var.id, sample=samp.id))
  seqClose(gds)
}


test_getGenotypeAlleles <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  alleles <- list(refChar(gds), altChar(gds, n=1), altChar(gds, n=2))
  for (i in 1:2) {
    for (j in 1:dim(geno)[2]) {
      for (a in 0:2) {
        x <- geno[i,j,] %in% a
        geno[i,j,x] <- alleles[[a+1]][x]
      }
    }
  }
  gc <- matrix(paste(geno[1,,], geno[2,,], sep="/"),
               nrow=dim(geno)[2], ncol=dim(geno)[3],
               dimnames=list(sample=NULL, variant=NULL))
  gc[gc == "NA/NA"] <- NA
  checkIdentical(gc, getGenotypeAlleles(gds, use.names=FALSE))
  seqClose(gds)
}

test_getGenotypeAlleles_phased <- function() {
  gds <- seqOpen(system.file("extdata", "gl_chr1.gds", package="SeqVarTools"))
  geno <- seqGetData(gds, "genotype")
  alleles <- list(refChar(gds), altChar(gds, n=1), altChar(gds, n=2))
  for (i in 1:2) {
    for (j in 1:dim(geno)[2]) {
      for (a in 0:2) {
        x <- geno[i,j,] %in% a
        geno[i,j,x] <- alleles[[a+1]][x]
      }
    }
  }
  gc <- matrix(paste(geno[1,,], geno[2,,], sep="|"),
               nrow=dim(geno)[2], ncol=dim(geno)[3],
               dimnames=list(sample=NULL, variant=NULL))
  gc[gc == "NA/NA"] <- NA
  checkIdentical(gc, getGenotypeAlleles(gds, use.names=FALSE))

  ## sort
  gc <- matrix(paste(pmin(geno[1,,], geno[2,,]),
                     pmax(geno[1,,], geno[2,,]), sep="/"),
               nrow=dim(geno)[2], ncol=dim(geno)[3],
               dimnames=list(sample=NULL, variant=NULL))
  gc[gc == "NA/NA"] <- NA
  checkIdentical(gc, getGenotypeAlleles(gds, use.names=FALSE, sort=TRUE))
  seqClose(gds)
}

test_getGenotypeAlleles_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  geno <- getGenotypeAlleles(gds)
  seqSetFilter(gds)
  checkIdentical(geno,
                 applyMethod(gds, getGenotypeAlleles, variant=var.id, sample=samp.id))
  seqClose(gds)
}

test_refDosage <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- getGenotype(gds)
  rd <- refDosage(gds)
  checkIdentical(geno %in% "0/0", rd %in% 2)
  checkIdentical(geno %in% c("0/1", "1/0"), rd %in% 1)
  checkIdentical(geno %in% "1/1", rd %in% 0)
  checkIdentical(is.na(geno), is.na(rd))
  seqClose(gds)
}

test_refDosage_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  geno <- refDosage(gds)
  seqSetFilter(gds)
  checkIdentical(geno,
                 applyMethod(gds, refDosage, variant=var.id, sample=samp.id))
  seqClose(gds)
}
