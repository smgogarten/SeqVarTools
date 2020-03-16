test_alleleFrequency_sum <- function() {
  gds <- SeqVarTools:::.testData()
  maxn <- max(nAlleles(gds))
  af <- matrix(nrow=SeqVarTools:::.nVar(gds), ncol=maxn)
  for (n in 1:maxn) af[,n] <- alleleFrequency(gds, n=(n-1))
  checkTrue(all(rowSums(af) == 1))
  seqClose(gds)
}

test_alleleFrequency_info <- function() {
  gds <- SeqVarTools:::.testData()
  ac <- seqGetData(gds, "annotation/info/AC")
  an <- seqGetData(gds, "annotation/info/AN")
  checkEquals(ac/an, alleleFrequency(gds, n=1))
  seqClose(gds)
}

test_alleleFrequency_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id, verbose=FALSE)
  af <- alleleFrequency(gds)
  seqSetFilter(gds, verbose=FALSE)
  checkIdentical(af,
                 applyMethod(gds, alleleFrequency, variant=var.id, sample=samp.id))
  seqClose(gds)
}

.testGdsXY <- function() {
    # make up file with sex chroms
    gds.fn <- tempfile()
    invisible(file.copy(seqExampleFileName("gds"), gds.fn))
    gds <- openfn.gds(gds.fn, readonly=FALSE)
    node <- index.gdsn(gds, "chromosome")
    compression.gdsn(node, "")
    chr <- read.gdsn(node)
    chr[chr == 1] <- "X"
    chr[chr == 2] <- "Y"
    write.gdsn(node, chr)
    closefn.gds(gds)
    seqOptimize(gds.fn, target="chromosome", verbose=FALSE)
    gds <- seqOpen(gds.fn)
    sample.id <- seqGetData(gds, "sample.id")
    set.seed(44); sex <- sample(c("M","F"), replace=TRUE, length(sample.id))
    df <- data.frame(sample.id, sex, stringsAsFactors=FALSE)
    SeqVarData(gds, sampleData=Biobase::AnnotatedDataFrame(df))
}

.cleanupGds <- function(gds) {
    fn <- seqSummary(gds, check="none", verbose=FALSE)$filename
    seqClose(gds)
    unlink(fn)
}


test_alleleFrequency_sex <- function() {
    svd <- .testGdsXY()
    sex <- sampleData(svd)$sex

    af <- alleleFrequency(svd)

    geno <- refDosage(svd, use.names=FALSE)
    chr <- chromWithPAR(svd)
    auto <- chr %in% 1:22
    checkEquals(0.5*colMeans(geno[,auto], na.rm=TRUE), af[auto])
    
    X <- chr == "X"
    female <- sex == "F"
    male <- sex == "M"
    F.count <- colSums(geno[female, X], na.rm=TRUE)
    F.nsamp <- colSums(!is.na(geno[female, X]))
    M.count <- 0.5*colSums(geno[male, X], na.rm=TRUE)
    M.nsamp <- colSums(!is.na(geno[male, X]))
    checkEquals((F.count + M.count)/(2*F.nsamp + M.nsamp), af[X])

    Y <- chr == "Y"
    checkEquals(0.5*colMeans(geno[male,Y], na.rm=TRUE), af[Y])

    # PAR
    checkTrue(all(chr[1:3] == "PAR"))
    checkEquals(0.5*colMeans(geno[,1:3], na.rm=TRUE), af[1:3])

    # names
    af <- alleleFrequency(svd, use.names=TRUE)
    checkEquals(as.character(seqGetData(svd, "variant.id")), names(af))
    
    .cleanupGds(svd)
}

test_alleleFrequency_nosex <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    # make sure there is no warning when we don't want to check sex
    options(warn=2)
    tmp <- alleleFrequency(gds, sex.adjust=FALSE)
    tmp <- alleleFrequency(gds, sex.adjust=TRUE) # no sex chrom in this data
    options(warn=1)
    seqClose(gds)
}


test_alleleCount_sex <- function() {
    svd <- .testGdsXY()
    sex <- sampleData(svd)$sex

    ac <- alleleCount(svd)

    geno <- refDosage(svd, use.names=FALSE)
    chr <- chromWithPAR(svd)
    auto <- chr %in% 1:22
    checkEquals(colSums(geno[,auto], na.rm=TRUE), ac[auto])
    
    X <- chr == "X"
    female <- sex == "F"
    male <- sex == "M"
    F.count <- colSums(geno[female, X], na.rm=TRUE)
    M.count <- 0.5*colSums(geno[male, X], na.rm=TRUE)
    checkEquals((F.count + M.count), ac[X])

    Y <- chr == "Y"
    checkEquals(0.5*colSums(geno[male,Y], na.rm=TRUE), ac[Y])

    # MAC
    mac <- minorAlleleCount(svd)
    ac.alt <- alleleCount(svd, n=1) + alleleCount(svd, n=2)
    minor <- ac < ac.alt
    checkEquals(ac[minor], mac[minor])
    checkEquals(ac.alt[!minor], mac[!minor])
    
    # PAR
    checkTrue(all(chr[1:3] == "PAR"))
    checkEquals(colSums(geno[,1:3], na.rm=TRUE), ac[1:3])

    # names
    ac <- alleleCount(svd, use.names=TRUE)
    checkEquals(as.character(seqGetData(svd, "variant.id")), names(ac))
    
    .cleanupGds(svd)
}


test_1KG_Y <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gds.fn <- system.file("extdata", "1KG_chrY.gds", package="SeqVarTools")
    
    gds <- seqOpen(gds.fn)
    sample.id <- seqGetData(gds, "sample.id")
    df <- data.frame(sample.id, sex="M", stringsAsFactors=FALSE)
    svd <- SeqVarData(gds, sampleData=Biobase::AnnotatedDataFrame(df))

    af <- alleleFrequency(svd)

    geno <- refDosage(svd, use.names=FALSE)
    chr <- chromWithPAR(gds)
    checkTrue(all(chr == "Y"))
    
    checkEquals(colMeans(geno, na.rm=TRUE), af)

    # AC
    ac <- alleleCount(svd)
    checkEquals(colSums(geno, na.rm=TRUE), ac)

    # MAC
    mac <- minorAlleleCount(svd)
    ac.alt <- alleleCount(svd, n=1) + alleleCount(svd, n=2) + alleleCount(svd, n=3) + alleleCount(svd, n=4) 
    minor <- ac < ac.alt
    checkEquals(ac[minor], mac[minor])
    checkEquals(ac.alt[!minor], mac[!minor])
    
    seqClose(gds)
}


test_alleleCount_allF <- function() {
    svd <- .testGdsXY()
    sex <- sampleData(svd)$sex
    seqSetFilter(svd, sample.sel=(sex == "F"), verbose=FALSE)

    ac <- alleleCount(svd)

    geno <- refDosage(svd, use.names=FALSE)
    chr <- chromWithPAR(svd)
    
    X <- chr == "X"
    F.count <- colSums(geno[, X], na.rm=TRUE)
    checkEquals((F.count), ac[X])

    Y <- chr == "Y"
    checkTrue(all(ac[Y] == 0))

    # MAC
    mac <- minorAlleleCount(svd)
    ac.alt <- alleleCount(svd, n=1) + alleleCount(svd, n=2)
    minor <- ac < ac.alt
    checkEquals(ac[minor], mac[minor])
    checkEquals(ac.alt[!minor], mac[!minor])
    
    .cleanupGds(svd)
}
