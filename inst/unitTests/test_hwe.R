
test_permute <- function() {
    x <- matrix(c(0,0,1,1,NA,NA), nrow=2)
    p <- SeqVarTools:::.permuteGenotypes(x)
    checkEquals(c(0,0,1,1), sort(as.vector(p[,1:2])))
    checkTrue(all(is.na(p[,3])))
}

test_count <- function() {
    gds <- seqOpen(seqExampleFileName("gds"))
    filt <- nAlleles(gds) == 2
    seqSetFilter(gds, variant.sel=filt)
    geno <- getGenotype(gds)
    counts <- data.frame(nAA=colSums(geno == "0/0", na.rm=TRUE),
                         nAa=colSums(geno == "0/1" | geno == "1/0", na.rm=TRUE),
                         naa=colSums(geno == "1/1", na.rm=TRUE),
                         row.names=1:ncol(geno))
    checkEquals(counts, SeqVarTools:::.countGenotypes(gds))
    perm <- SeqVarTools:::.countGenotypes(gds, permute=TRUE)
    checkEquals(rowSums(counts), rowSums(perm), checkNames=FALSE)
    checkEquals(2*counts$nAA + counts$nAa, 2*perm$nAA + perm$nAa, checkNames=FALSE)
    seqClose(gds)
}

test_hwe <- function() {
    gds <- seqOpen(seqExampleFileName("gds"))
    af <- alleleFrequency(gds)
    biallelic <- nAlleles(gds) == 2
    mono <- af %in% c(0,1)
    hw <- hwe(gds, permute=FALSE)
    checkEquals(mono | !biallelic, is.na(hw$p), checkNames=FALSE)
    
    filt <- biallelic & !mono
    seqSetFilter(gds, variant.sel=filt)
    hw <- hwe(gds, permute=FALSE)
    checkEquals(af[filt], hw$afreq)
    checkEquals(inbreedCoeff(gds), hw$f)

    hw <- hwe(gds, permute=TRUE)
    checkEquals(af[filt], hw$afreq)
    seqClose(gds)
}
