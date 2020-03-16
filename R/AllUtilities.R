
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

.rangesToSel <- function(gdsobj, ranges) {
  gds.ranges <- granges(gdsobj)
  # don't warn if no sequence levels in common
  queryHits(suppressWarnings(findOverlaps(gds.ranges, ranges)))
}

.ploidy <- function(gdsobj) {
  seqSummary(gdsobj, "genotype", check="none", verbose=FALSE)$dim[1L]
}

.nSamp <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$sample.sel)
  seqSummary(gdsobj, "genotype", check="none", verbose=FALSE)$seldim[2L]
}

.nVar <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$variant.sel)
  seqSummary(gdsobj, "genotype", check="none", verbose=FALSE)$seldim[3L]
}

.nSampUnfiltered <- function(gdsobj) {
  seqSummary(gdsobj, "sample.id", check="none", verbose=FALSE)
}

.nVarUnfiltered <- function(gdsobj) {
  seqSummary(gdsobj, "variant.id", check="none", verbose=FALSE)
}

.nSampObserved <- function(gdsobj) {
  ns <- .nSamp(gdsobj)
  if (ns > 0) {
      return(ns * (1-seqMissing(gdsobj)))
  } else {
      return(rep(0, .nVar(gdsobj)))
  }
}

.emptySampFilter <- function(x, verbose=FALSE) {
    seqSetFilter(x, sample.sel=raw(.nSampUnfiltered(x)), action="push+set", verbose=verbose)
}

.emptyVarFilter <- function(x, verbose=FALSE) {
    seqSetFilter(x, variant.sel=raw(.nVarUnfiltered(x)), action="push+set", verbose=verbose)
}

.emptyDim <- function(x) {
    .nSamp(x) == 0 | .nVar(x) == 0
}

.emptyGenoMatrix <- function(x, use.names=FALSE) {
    m <- matrix(nrow=.nSamp(x), ncol=.nVar(x),
                dimnames=list(sample=NULL, variant=NULL))
    if (use.names) .applyNames(x, m) else m
}

.applyNames <- function(gdsobj, var) {
  if ("sample" %in% names(dimnames(var)))
    dimnames(var)$sample <- seqGetData(gdsobj, "sample.id")
  if ("variant" %in% names(dimnames(var)))
    dimnames(var)$variant <- seqGetData(gdsobj, "variant.id")
  var
}

.parseRefAllele <- function(x) {
  endRef <- regexpr(",", x, fixed=TRUE) - 1L
  noAlt <- endRef < 0L
  if (any(noAlt))
      endRef[noAlt] <- nchar(x[noAlt])
  substr(x, 1L, endRef)
}

.parseAltAllele <- function(x, n=0) {
  if (n == 0) {
    ## if n=0, return string with all ALT alleles
    begAlt <- regexpr(",", x, fixed=TRUE) + 1L
    noAlt <- begAlt == 0L
    if (any(noAlt))
        begAlt[noAlt] <- nchar(x[noAlt]) + 1L
    substr(x, begAlt, nchar(x))
  } else {
    ## if n>0, return nth ALT allele
    unlist(lapply(strsplit(x, ",", fixed=TRUE), function(x) x[n+1]),
           use.names=FALSE)
  }
}

## .parseNumAlleles <- function(x) {
##   #unlist(lapply(strsplit(x, ",", fixed=TRUE), length), use.names=FALSE)
##   str_count(x, ",") + 1L
## }

.maxAlleleLength <- function(x) {
  a <- gregexpr("[ACGT]+", x)
  unlist(lapply(a, function(y) max(attr(y, "match.length"))), use.names=FALSE)
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

.isTransition <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("C","T")) |
  (ref %in% c("A","G") & alt %in% c("A","G"))
}

.isTransversion <- function(ref, alt) {
  (ref %in% c("C","T") & alt %in% c("A","G")) |
  (ref %in% c("A","G") & alt %in% c("C","T"))
}

## number of homozygote reference genotypes
.nHomRef <- function(x) {
    sum(x[1,] == 0 & x[2,] == 0, na.rm=TRUE)
}

## number of heterozyotes with one reference allele
.nHetRef <- function(x) {
    sum((x[1,] == 0 & x[2,] != 0) |
        (x[1,] != 0 & x[2,] == 0), na.rm=TRUE)
}

## number of genotypes with no reference alleles
.nHomAlt <- function(x) {
    sum(x[1,] != 0 & x[2,] != 0, na.rm=TRUE)
}

.permuteGenotypes <- function(x) {
    ## get subset of matrix with no missing values
    ind <- colSums(is.na(x)) == 0
    nm <- x[,ind,drop=FALSE]
    ## permute alleles
    nm <- matrix(sample(nm), nrow=nrow(nm), ncol=ncol(nm))
    ## replace non-missing genotypes
    x[,ind] <- nm
    x
}

.countGenotypes <- function(gdsobj, permute=FALSE) {
    n <- seqApply(gdsobj, "genotype", function(x) {
        if (permute) x <- .permuteGenotypes(x)
        c(nAA=.nHomRef(x), nAa=.nHetRef(x), naa=.nHomAlt(x))
    }, margin="by.variant", as.is="list")
    as.data.frame(do.call(rbind, n))
}

.subjectByQuery <- function(query, subject, hits.only=FALSE) {
    ol <- findOverlaps(query, subject)
    if (hits.only) {
        hits <- unique(queryHits(ol))
    } else {
        hits <- seq_along(query)
    }
    subj.hits <- lapply(hits, function(h) {
        subjectHits(ol)[queryHits(ol) == h]
    })
    list(queryHits=hits, subjectHits=subj.hits)
}

# behaves like subsetByOverlaps, except removes query ranges that have
# duplicate sets of overlapping subject ranges
.subsetByUniqueOverlaps <- function(query, subject) {
    # look for any subject ranges that have identical query ranges
    # we want to eliminate these
    sbq <- .subjectByQuery(query, subject, hits.only=TRUE)
    query[sbq$queryHits[!duplicated(sbq$subjectHits)]]
}

.testData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- seqExampleFileName("gds")
    seqOpen(gdsfile)
}

.testSeqVarData <- function() {
    SeqVarData(.testData())
}

.testDosageData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "gl_chr1.gds", package="SeqVarTools")
    seqOpen(gdsfile)
}

.alleleCount <- function(gdsobj, n=0) {
    freq <- seqAlleleFreq(gdsobj, ref.allele=n)
    nsamp <- .nSamp(gdsobj)*(1-seqMissing(gdsobj))
    2*freq*nsamp
}
