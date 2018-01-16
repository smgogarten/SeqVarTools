
test_altDosage_sparse <- function() {
  gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
  gds <- seqOpen(gdsfile)
  d1 <- altDosage(gds, sparse=FALSE)
  d2 <- altDosage(gds, sparse=TRUE)
  checkTrue(is.matrix(d1))
  checkTrue(is(d2, "dgCMatrix"))
  checkTrue(all(is.na(d1) == is.na(d2)))
  checkTrue(all(is.na(d2) | d1 == d2))
  checkEquals(dimnames(d1), dimnames(d2))
  seqClose(gds)
}

test_expandedAltDosage_sparse <- function() {
  gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
  gds <- seqOpen(gdsfile)
  d1 <- expandedAltDosage(gds, sparse=FALSE)
  d2 <- expandedAltDosage(gds, sparse=TRUE)
  checkTrue(is.matrix(d1))
  checkTrue(is(d2, "dgCMatrix"))
  checkTrue(all(is.na(d1) == is.na(d2)))
  checkTrue(all(is.na(d2) | d1 == d2))
  checkEquals(dimnames(d1), dimnames(d2))
  seqClose(gds)
}
