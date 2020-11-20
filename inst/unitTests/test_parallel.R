test_parallel_bysample  <- function() {
    gds <- SeqVarTools:::.testData()
    
    ic1 <- inbreedCoeff(gds, margin="by.sample", parallel=FALSE)
    ic2 <- inbreedCoeff(gds, margin="by.sample", parallel=TRUE)
    checkEquals(ic1, ic2)

    miss1 <- missingGenotypeRate(gds, margin="by.sample", parallel=FALSE)
    miss2 <- missingGenotypeRate(gds, margin="by.sample", parallel=TRUE)
    checkEquals(miss1, miss2)
    
    het1 <- heterozygosity(gds, margin="by.sample", parallel=FALSE)
    het2 <- heterozygosity(gds, margin="by.sample", parallel=TRUE)
    checkEquals(het1, het2)

    seqClose(gds)
}

test_parallel_byvariant  <- function() {
    gds <- SeqVarTools:::.testData()
    
    ic1 <- inbreedCoeff(gds, margin="by.variant", parallel=FALSE)
    ic2 <- inbreedCoeff(gds, margin="by.variant", parallel=TRUE)
    checkEquals(ic1, ic2)

    miss1 <- missingGenotypeRate(gds, margin="by.variant", parallel=FALSE)
    miss2 <- missingGenotypeRate(gds, margin="by.variant", parallel=TRUE)
    checkEquals(miss1, miss2)

    hwe1 <- hwe(gds, parallel=FALSE)
    hwe2 <- hwe(gds, parallel=TRUE)
    checkEquals(hwe1, hwe2)

    het1 <- heterozygosity(gds, margin="by.variant", parallel=FALSE)
    het2 <- heterozygosity(gds, margin="by.variant", parallel=TRUE)
    checkEquals(het1, het2)
    
    seqClose(gds)
}
