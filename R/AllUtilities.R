
.applyMethod <- function(gdsobj, FUN, variant.id=NULL, sample.id=NULL, ...) {
  seqSetFilter(gdsobj, sample.id=sample.id, variant.id=variant.id, action="push")
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
    unlist(lapply(strsplit(x, ",", fixed=TRUE), function(x) x[n+1]))
  }
}

.parseNumAlleles <- function(x) {
  unlist(lapply(strsplit(x, ",", fixed=TRUE), length))
}

.maxAlleleLength <- function(x) {
  a <- gregexpr("[ACGT]+", x)
  unlist(lapply(a, function(y) max(attr(y, "match.length"))))
}
  
.isTransition <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("C","T")) |
  (ref %in% c("A","G") & alt %in% c("A","G"))
}

.isTransversion <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("A","G")) |
  (ref %in% c("A","G") & alt %in% c("C","T"))
}

.parseVariableLength <- function(x) {
  if (all(x$length == 1)) {
    x$data
  } else {
    var <- array(dim=c(max(x$length), nrow(x$data), length(x$length)))

    ## assign each element of length to an index of first array dimension
    n.ind <- rep(NA, ncol(x$data))
    j <- 1
    for (i in 1:length(x$length)) {
      len <- x$length[i]
      if (len > 0) {
        n.ind[j:(j + len - 1)] <- 1:len
        j <- j + len
      }
    }

    ## for each index of first array dimension, get values
    for (n in 1:dim(var)[1]) {
      var.ind <- which(x$length >= n)
      var[n,,var.ind] <- x$data[,n.ind == n]
    }

    ## if first array dimension is 1, simplify to a matrix
    if (dim(var)[1] == 1) {
      var <- var[1,,]
    }
    var
  }
}
