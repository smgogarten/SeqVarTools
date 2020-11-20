
setMethod("inbreedCoeff",
          "SeqVarGDSClass",
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE, parallel=FALSE) {           
            margin <- match.arg(margin)
            if (margin == "by.variant") {
              counts <- .countGenotypes(gdsobj, parallel=parallel)
              f <- .f(counts)
              if (use.names)
                names(f) <- seqGetData(gdsobj, "variant.id")
              f
            } else {
              nsamp <- .nSamp(gdsobj)

              ## need global variables
              n <- integer(nsamp)
              s <- double(nsamp)

              ## apply the function variant by variant
              seqApply(gdsobj, "genotype", function(x) {
                p <- mean(x==0L, na.rm=TRUE)     # allele frequency
                g <- colSums(x==0L)              # genotypes: # of reference allele
                d <- (g*g - g*(1 + 2*p) + 2*p*p) / (2*p*(1-p))
                n <<- n + is.finite(d)           # output to the global variable "n"
                d[!is.finite(d)] <- 0
                s <<- s + d                      # output to the global variable "s"
              }, margin="by.variant", as.is="none")

              ## output
              ic <- s / n
              if (use.names)
                names(ic) <- seqGetData(gdsobj, "sample.id")
              ic
            }
          })
