test_refFrac_excep <- function() {
    ## no AD field
    gds <- seqOpen(seqExampleFileName("gds"))
    checkException(refFrac(gds))
    seqClose(gds)
}

test_refFrac <- function() {
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.sel=isSNV(gds))
    ad <- getVariableLengthData(gds, "annotation/format/AD")
    rf <- ad[1,,] / (ad[1,,] + ad[2,,])
    checkEquals(rf, refFrac(gds))

    geno <- seqGetData(gds, "genotype")
    het <- geno[1,,] != geno[2,,]
    rf[!het | is.na(het)] <- NA
    checkEquals(colMeans(rf, na.rm=TRUE), refFracOverHets(gds), checkNanes=FALSE)
    seqClose(gds)
}
