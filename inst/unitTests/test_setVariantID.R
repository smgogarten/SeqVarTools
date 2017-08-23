test_setVariantID <- function() {
  gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
  oldfile <- system.file("extdata", "gl_chr1.gds", package="SeqVarTools")
  newfile <- tempfile()
  file.copy(oldfile, newfile)

  # variant length
  gds <- seqOpen(newfile)
  var.len <- length(seqGetData(gds, "variant.id"))
  seqClose(gds)

  # integer id
  int.id <- seq(101, 100+var.len)
  setVariantID(newfile, int.id)
  gds <- seqOpen(newfile)
  checkIdentical(int.id, seqGetData(gds, "variant.id"))
  seqClose(gds)
 
  # character id
  char.id <- letters[1:var.len]
  setVariantID(newfile, char.id)
  gds <- seqOpen(newfile)
  checkIdentical(char.id, seqGetData(gds, "variant.id"))
  seqClose(gds)

  # non-unique
  rep.id <- rep("a", var.len)
  checkException(setVariantID(newfile, rep.id))

  # wrong length
  long.id <- 1:(var.len+1)
  checkException(setVariantID(newfile, long.id))
 
  unlink(newfile)
}
