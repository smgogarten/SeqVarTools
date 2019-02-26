test_applyNames <- function() {
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
  geno <- seqGetData(gds, "genotype")
  var <- (!is.na(geno[1,,]) & geno[1,,] > 0) | (!is.na(geno[2,,]) & geno[2,,] > 0)
  dimnames(var) <- list(sample=NULL, variant=NULL)
  checkIdentical(var, isVariant(gds))
  seqClose(gds)
}

test_isVariant_apply <- function() {
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
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
  gds <- SeqVarTools:::.testData()
  geno <- getGenotype(gds)
  rd <- refDosage(gds)
  checkIdentical(geno %in% "0/0", rd %in% 2)
  checkIdentical(geno %in% c("0/1", "1/0"), rd %in% 1)
  checkIdentical(geno %in% "1/1", rd %in% 0)
  checkIdentical(is.na(geno), is.na(rd))
  seqClose(gds)
}

test_refDosage_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  geno <- refDosage(gds)
  seqSetFilter(gds)
  checkIdentical(geno,
                 applyMethod(gds, refDosage, variant=var.id, sample=samp.id))
  seqClose(gds)
}

test_altDosage <- function() {
  gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
  gds <- seqOpen(gdsfile)
  geno <- getGenotype(gds)
  d <- altDosage(gds)
  hasref <- geno %in% paste0("0/", 0:(max(nAlleles(gds))-1))
  hasref[is.na(geno)] <- NA
  hasalt <- !(geno %in% "0/0")
  hasalt[is.na(geno)] <- NA
  checkEquals(!hasref, as.vector(d == 2))
  checkEquals(hasalt & hasref, as.vector(d == 1))
  checkEquals(geno == "0/0", d == 0)
  checkIdentical(is.na(geno), is.na(d))
  checkIdentical(refDosage(gds), 2-d)
  seqClose(gds)
}

## former alleleDosage method for list, summing over all requested alleles
.alleleDosageSum <- function(gdsobj, n, use.names=TRUE) {
    stopifnot(length(n) == SeqVarTools:::.nVar(gdsobj))
    d <- seqApply(gdsobj, "genotype",
                  function(index, x) {
                      tmp <- matrix(x %in% n[[index]], ncol=ncol(x),
                                    nrow=nrow(x))
                      tmp[is.na(x)] <- NA
                      colSums(tmp)
                  },
                  margin="by.variant", as.is="list",
                  var.index="relative")
    d <- matrix(unlist(d, use.names=FALSE), ncol=length(d),
                dimnames=list(sample=NULL, variant=NULL))
    if (use.names) SeqVarTools:::.applyNames(gdsobj, d) else d
}

test_alleleDosage <- function() {
  gds <- SeqVarTools:::.testData()
  checkEquals(refDosage(gds), alleleDosage(gds, n=0))
  checkEquals(altDosage(gds), alleleDosage(gds, n=1))
  seqClose(gds)
  
  gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
  gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
  gds <- seqOpen(gdsfile)
  checkEquals(refDosage(gds), alleleDosage(gds, n=0))
  checkEquals(altDosage(gds), 2-alleleDosage(gds, n=0))

  ## different alleles for each variant
  nalleles <- nAlleles(gds)
  n <- sapply(nalleles, function(x) sample(0:(x-1), 1))
  geno <- getGenotype(gds)
  cnt <- matrix(nrow=nrow(geno), ncol=ncol(geno), dimnames=dimnames(geno))
  for (i in 1:nrow(geno)) {
      cnt[i,] <- stringr::str_count(geno[i,], as.character(n))
  }
  checkEquals(cnt, alleleDosage(gds, n))

  ## n is a list
  n <- lapply(nalleles, function(x) sample(0:(x-1), sample(1:x, 1)))
  cnt <- list()
  for (i in 1:ncol(geno)) {
  ##     tmp <- lapply(n[[i]], function(x) stringr::str_count(geno[,i], as.character(x)))
  ##     cnt[,i] <- colSums(do.call(rbind, tmp))
      cnt[[i]] <- do.call(cbind, lapply(n[[i]], function(x) stringr::str_count(geno[,i], as.character(x))))
      dimnames(cnt[[i]]) <- list(sample=rownames(geno), allele=n[[i]])
  }
  ad <- alleleDosage(gds, n)
  checkEquals(cnt, ad, checkNames=FALSE)

  adsum <- do.call(cbind, lapply(ad, function(x) rowSums(x)))
  names(dimnames(adsum)) <- c("sample", "variant")
  checkEquals(.alleleDosageSum(gds, n), adsum)
  
  ## invalid allele
  checkException(alleleDosage(gds, n=10))
  checkException(alleleDosage(gds, n=c(1,0)))
  
  seqClose(gds)
}

test_expandedAltDosage <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    ead <- expandedAltDosage(gds)
    checkEquals(sum(nAlleles(gds) - 1), ncol(ead)) 

    seqSetFilter(gds, variant.sel=(nAlleles(gds) == 3), verbose=FALSE)
    ad <- lapply(1:2, function(n) alleleDosage(gds, n=n))
    ead <- expandedAltDosage(gds)
    checkEquals(ad[[1]], ead[,c(TRUE,FALSE)])
    checkEquals(ad[[2]], ead[,c(FALSE,TRUE)])

    seqResetFilter(gds, verbose=FALSE)
    seqSetFilter(gds, variant.sel=(nAlleles(gds) == 4), verbose=FALSE)
    ad <- lapply(1:3, function(n) alleleDosage(gds, n=n))
    ead <- expandedAltDosage(gds)
    checkEquals(ad[[1]], ead[,c(TRUE,FALSE,FALSE)])
    checkEquals(ad[[2]], ead[,c(FALSE,TRUE,FALSE)])
    checkEquals(ad[[3]], ead[,c(FALSE,FALSE,TRUE)])

    seqClose(gds)
}

test_expandedVariantIndex <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    ind <- expandedVariantIndex(gds)
    checkEquals(sum(nAlleles(gds) - 1), length(ind)) 
    ead <- expandedAltDosage(gds)
    checkEquals(ncol(ead), length(ind))
    checkEquals(colnames(ead), as.character(seqGetData(gds, "variant.id")[ind]))
    checkEquals(as.vector(table(nAlleles(gds))), as.vector(table(table(ind))))
    seqClose(gds)
}

test_expandedVariantIndex_noAlt <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "gl_chr1.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    ind <- expandedVariantIndex(gds)
    checkEquals(1:9, ind) 
    seqClose(gds)
}

test_variantInfo <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    
    v <- variantInfo(gds, alleles=FALSE, expanded=FALSE)
    checkEquals(SeqVarTools:::.nVar(gds), nrow(v))
    
    v <- variantInfo(gds, alleles=TRUE, expanded=FALSE)
    checkTrue(all(c("ref", "alt") %in% names(v)))
    
    v <- variantInfo(gds, alleles=FALSE, expanded=TRUE)
    na <- table(nAlleles(gds))
    ai <- sapply(1:length(na), function(i) sum(na[i:length(na)]))
    checkEquals(ai, as.vector(table(v$allele.index)))

    v <- variantInfo(gds, alleles=TRUE, expanded=TRUE)
    checkEquals(ai, as.vector(table(v$allele.index)))
    for (i in 1:6) {
        ac <- altChar(gds, n=i)
        checkEquals(ac[!is.na(ac)], v$alt[v$allele.index == i])
    }
    
    seqClose(gds)
}

test_empty <- function() {
    gds <- SeqVarTools:::.testData()
    SeqVarTools:::.emptyVarFilter(gds)
    checkEquals(0, ncol(getGenotype(gds)))
    checkEquals(0, ncol(getGenotypeAlleles(gds)))
    checkEquals(0, ncol(refDosage(gds)))
    checkEquals(0, ncol(altDosage(gds)))
    checkEquals(0, ncol(expandedAltDosage(gds)))
    checkEquals(0, ncol(alleleDosage(gds, n=0)))
    checkEquals(0, length(expandedVariantIndex(gds)))
    checkEquals(0, nrow(variantInfo(gds)))
    checkEquals(0, nrow(variantInfo(gds, expanded=TRUE)))
    checkEquals(0, nrow(variantInfo(gds, alleles=FALSE, expanded=TRUE)))
    seqClose(gds)
}
