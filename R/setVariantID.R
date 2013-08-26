
setVariantID <- function(gdsfile, variant.id) {
  stopifnot(file.exists(gdsfile))
  if (length(unique(variant.id)) != length(variant.id)) {
    stop("variant.id must be unique")
  }

  gdsobj <- openfn.gds(gdsfile, readonly=FALSE)
  node <- index.gdsn(gdsobj, "variant.id")
  curr.id <- read.gdsn(node)
  if (length(curr.id) != length(variant.id)) {
    closefn.gds(gdsobj)
    stop(paste("variant.id must have length", length(curr.id)))
  }

  ## delete existing node and create a new one,
  ## as data type may be different
  compress <- objdesp.gdsn(node)$compress
  delete.gdsn(node)
  add.gdsn(gdsobj, "variant.id", variant.id,
           compress=compress, closezip=TRUE)
  closefn.gds(gdsobj)
}
