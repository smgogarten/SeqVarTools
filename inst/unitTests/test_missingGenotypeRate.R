test_missingGenotypeRate <- function() {
  gds <- SeqVarTools:::.testData()
  geno <- seqGetData(gds, "genotype")
  checkIdentical(colSums(is.na(geno[1,,])) / dim(geno)[2],
                 missingGenotypeRate(gds, "by.variant"))
  checkIdentical(rowSums(is.na(geno[1,,])) / dim(geno)[3],
                 missingGenotypeRate(gds, "by.sample"))
  seqClose(gds)
}

test_missingGenotypeRate_apply <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  miss.var <- missingGenotypeRate(gds, margin="by.variant")
  miss.samp <- missingGenotypeRate(gds, margin="by.sample")
  seqSetFilter(gds)
  checkIdentical(miss.var,
                 applyMethod(gds, missingGenotypeRate, variant=var.id, sample=samp.id,
                             margin="by.variant"))
  checkIdentical(miss.samp,
                 applyMethod(gds, missingGenotypeRate, variant=var.id, sample=samp.id,
                             margin="by.sample"))
  seqClose(gds)
}
