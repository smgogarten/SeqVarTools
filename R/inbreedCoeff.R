
setMethod("inbreedCoeff",
          "SeqVarGDSClass",
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE) {           
            margin <- match.arg(margin)
            if (margin == "by.variant") {
              nhet <- seqApply(gdsobj, "genotype",
                               function(x) {sum((x[1,] == 0 & x[2,] == 1) |
                                                (x[1,] == 1 & x[2,] == 0), na.rm=TRUE)},
                               margin="by.variant", as.is="integer")
              ntot <- seqApply(gdsobj, "genotype",
                               function(x) {sum(!is.na(x[1,]))},
                               margin="by.variant", as.is="integer")
              afreq <- seqApply(gdsobj, "genotype",
                                function(x) {mean(x == 0, na.rm=TRUE)},
                               margin="by.variant", as.is="double")
              exp.het <- 2*afreq*(1-afreq)*ntot
              f <- 1 - (nhet/exp.het)
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
                p <- mean(x==0, na.rm=TRUE)      # allele frequency
                g <- colSums(x==0)               # genotypes: # of reference allele
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
