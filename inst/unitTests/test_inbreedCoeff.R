
test_inbreedCoeff_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  ic.var <- inbreedCoeff(gds, margin="by.variant")
  ic.samp <- inbreedCoeff(gds, margin="by.sample")
  seqSetFilter(gds)
  checkIdentical(ic.var,
                 applyMethod(gds, inbreedCoeff, variant=var.id, sample=samp.id,
                             margin="by.variant"))
  checkIdentical(ic.samp,
                 applyMethod(gds, inbreedCoeff, variant=var.id, sample=samp.id,
                             margin="by.sample"))
  seqClose(gds)
}

