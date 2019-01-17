test_dosage <- function() {
    gds <- SeqVarTools:::.testDosageData()
    dos <- imputedDosage(gds)
    checkTrue(is(dos, "matrix"))
    checkEquals(rownames(dos), seqGetData(gds, "sample.id"))
    checkEquals(colnames(dos), as.character(seqGetData(gds, "variant.id")))
    checkTrue(all(dos >= 0 & dos <= 2))
    seqClose(gds)
}
