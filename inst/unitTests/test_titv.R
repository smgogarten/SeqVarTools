test_titv_bysamp <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  ref <- refChar(gds)
  alt <- altChar(gds)
  ti <- SeqVarTools:::.isTransition(ref, alt)
  tv <- SeqVarTools:::.isTransversion(ref, alt)
  var <- isVariant(gds)
  timat <- matrix(ti, nrow=nrow(var), ncol=ncol(var), byrow=TRUE)
  tvmat <- matrix(tv, nrow=nrow(var), ncol=ncol(var), byrow=TRUE)
  checkIdentical(rowSums(timat & var) / rowSums(tvmat & var),
                 titv(gds, by.sample=TRUE))
  seqClose(gds)
}

test_titv_bysamp_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  ti.tv <- titv(gds, by.sample=TRUE)
  seqSetFilter(gds)
  checkIdentical(ti.tv,
                 applyMethod(gds, titv, variant=var.id, sample=samp.id,
                             by.sample=TRUE))
  seqClose(gds)
}

test_titv_byvar_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  ti.tv <- titv(gds)
  seqSetFilter(gds)
  checkIdentical(ti.tv,
                 applyMethod(gds, titv, variant=var.id, sample=samp.id))
  seqClose(gds)
}
