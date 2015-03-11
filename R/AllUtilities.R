
.applyMethod <- function(gdsobj, FUN, variant.id=NULL, sample.id=NULL, ...) {
  seqSetFilter(gdsobj, sample.id=sample.id, variant.id=variant.id, action="push+set")
  result <- FUN(gdsobj, ...)
  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
  result
}

.rangesToID <- function(gdsobj, ranges) {
  gds.ranges <- granges(gdsobj, id=seqGetData(gdsobj, "variant.id"))
  mcols(subsetByOverlaps(gds.ranges, ranges))$id
}

.nSamp <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$sample.sel)
  seqSummary(gdsobj, "genotype")$seldim[1]
}

.nVar <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$variant.sel)
  seqSummary(gdsobj, "genotype")$seldim[2]
}

.applyNames <- function(gdsobj, var) {
  if ("sample" %in% names(dimnames(var)))
    dimnames(var)$sample <- seqGetData(gdsobj, "sample.id")
  if ("variant" %in% names(dimnames(var)))
    dimnames(var)$variant <- seqGetData(gdsobj, "variant.id")
  var
}

.parseRefAllele <- function(x) {
  SeqArray:::.refAllele(x)
}

.parseAltAllele <- function(x, n=0) {
  if (n == 0) {
    ## if n=0, return string with all ALT alleles
    SeqArray:::.altAllele(x)
  } else {
    ## if n>0, return nth ALT allele
    unlist(lapply(strsplit(x, ",", fixed=TRUE), function(x) x[n+1]),
           use.names=FALSE)
  }
}

.parseNumAlleles <- function(x) {
  unlist(lapply(strsplit(x, ",", fixed=TRUE), length), use.names=FALSE)
}

.maxAlleleLength <- function(x) {
  a <- gregexpr("[ACGT]+", x)
  unlist(lapply(a, function(y) max(attr(y, "match.length"))), use.names=FALSE)
}
  
.isTransition <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("C","T")) |
  (ref %in% c("A","G") & alt %in% c("A","G"))
}

.isTransversion <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("A","G")) |
  (ref %in% c("A","G") & alt %in% c("C","T"))
}
