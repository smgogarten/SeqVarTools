test_applyMethod_util <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- seqGetData(gds, "variant.id") 
  checkIdentical(var.id,
                 SeqVarTools:::.applyMethod(gds, seqGetData, var.name="variant.id"))
  checkIdentical(var.id[1:5],
                 SeqVarTools:::.applyMethod(gds, seqGetData, variant=var.id[1:5], var.name="variant.id")) 
  samp.id <- seqGetData(gds, "sample.id")
  checkIdentical(samp.id,
                 SeqVarTools:::.applyMethod(gds, seqGetData, var.name="sample.id"))
  checkIdentical(samp.id[1:5],
                 SeqVarTools:::.applyMethod(gds, seqGetData, sample=samp.id[1:5], var.name="sample.id"))
  seqClose(gds)
}

test_rangesToID <- function() {
  gds <- SeqVarTools:::.testData()
  chrom <- seqGetData(gds, "chromosome")
  var.id <- seqGetData(gds, "variant.id")[chrom == 1]
  pos <- seqGetData(gds, "position")[chrom == 1]
  ranges <- GRanges(seqnames="1", IRanges(start=min(pos), end=max(pos)))
  checkIdentical(var.id, SeqVarTools:::.rangesToID(gds, ranges))
  seqClose(gds)
}
  
test_applyMethod_vector <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- seqGetData(gds, "variant.id") 
  checkIdentical(var.id[1:5],
                 applyMethod(gds, seqGetData, variant=var.id[1:5], var.name="variant.id"))
  seqClose(gds)
}

test_applyMethod_ranges <- function() {
  gds <- SeqVarTools:::.testData()
  chrom <- seqGetData(gds, "chromosome")
  var.id <- seqGetData(gds, "variant.id")[chrom == 1]
  pos <- seqGetData(gds, "position")[chrom == 1]
  ranges <- GRanges(seqnames="1", IRanges(start=min(pos), end=max(pos)))
  checkIdentical(var.id,
                 applyMethod(gds, seqGetData, variant=ranges, var.name="variant.id"))
  seqClose(gds)
}

test_applyMethod_missing <- function() {
  gds <- SeqVarTools:::.testData()
  var.id <- seqGetData(gds, "variant.id") 
  checkIdentical(var.id,
                 applyMethod(gds, seqGetData, var.name="variant.id"))
  seqClose(gds)
}

test_applyMethod_pushpop <- function() {
  gds <- SeqVarTools:::.testData()
  sample.id <- seqGetData(gds, "sample.id")
  seqSetFilter(gds, sample.id=sample.id[1:5], variant.id=1:10)
  samp.orig <- seqGetData(gds, "sample.id")
  var.orig <- seqGetData(gds, "variant.id")
  applyMethod(gds, print, sample=sample.id[6:10], variant=11:20)
  checkIdentical(samp.orig, seqGetData(gds, "sample.id"))
  checkIdentical(var.orig, seqGetData(gds, "variant.id"))
  seqClose(gds)
}
