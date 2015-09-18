.nHomRef <- function(x) {
    sum(x[1,] == 0 & x[2,] == 0, na.rm=TRUE)
}

.nHet <- function(x) {
    sum((x[1,] == 0 & x[2,] != 0) |
        (x[1,] != 0 & x[2,] == 0), na.rm=TRUE)
}

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
        c(nAA=.nHomRef(x), nAa=.nHet(x), naa=.nHomAlt(x))
    }, margin="by.variant", as.is="list")
    counts <- as.data.frame(do.call(rbind, n))
}

.f <- function(counts) {
    nhet <- counts$nAa
    ntot <- rowSums(counts)
    afreq <- (2*counts$nAA + counts$nAa) / (2*ntot)
    exp.het <- 2*afreq*(1-afreq)*ntot
    1 - (nhet/exp.het)
}

setMethod("hwe",
          "SeqVarGDSClass",
          function(gdsobj, permute=FALSE) {
              counts <- .countGenotypes(gdsobj, permute)

              afreq <- (2*counts$nAA + counts$nAa) / (2*rowSums(counts))
              p <- HWExact(counts)
              f <- .f(counts)

              ## set non-biallelic or monomorphic variants to NA
              sel <- nAlleles(gdsobj) != 2 |
                  counts$nAa + counts$naa == 0 |
                      counts$nAa + counts$nAA == 0
              p[sel] <- NA
              f[sel] <- NA

              variant.id <- seqGetData(gdsobj, "variant.id")
              cbind(variant.id, counts, afreq, p, f)
          })
