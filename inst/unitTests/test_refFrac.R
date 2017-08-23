test_refFrac_excep <- function() {
    ## no AD field
    gds <- seqOpen(seqExampleFileName("gds"))
    checkException(refFrac(gds))
    seqClose(gds)
}

test_refFrac <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.sel=isSNV(gds))
    ad <- getVariableLengthData(gds, "annotation/format/AD")
    rf <- ad[1,,] / (ad[1,,] + ad[2,,])
    checkEquals(rf, refFrac(gds))

    geno <- seqGetData(gds, "genotype")
    het <- geno[1,,] != geno[2,,]
    rf[!het | is.na(het)] <- NA
    checkEquals(colMeans(rf, na.rm=TRUE), refFracOverHets(gds), checkNames=FALSE)

    ## compare with allelic balance from VCF
    ab <- getVariableLengthData(gds, "annotation/format/AB")
    rf[is.na(ab)] <- NA
    checkEquals(rf, ab, tolerance=0.01)
    
    seqClose(gds)
}
