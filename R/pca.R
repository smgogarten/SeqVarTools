
setMethod("pca",
          "SeqVarGDSClass",
          function(gdsobj, eigen.cnt=32) {
            nsamp <- .nSamp(gdsobj)
            if (eigen.cnt > nsamp) eigen.cnt <- nsamp

            ## genetic correlation matrix
            genmat <- local({

              ## need a global variable (only available in the bracket of "local")
              s <- matrix(0.0, nrow=nsamp, ncol=nsamp)

              ## apply the function variant by variant
              seqApply(gdsobj, "genotype", function(x) {
                g <- (x==0L)                  # indicator of reference allele
                p <- mean(g, na.rm=TRUE)      # allele frequency
                g2 <- colSums(g) - 2*p        # genotypes adjusted by allele frequency
                m <- (g2 %o% g2) / (p*(1-p))  # genetic correlation matrix
                m[!is.finite(m)] <- 0         # correct missing values
                s <<- s + m                   # output to the global variable "s"
              }, margin="by.variant", as.is="none")

              ## output, scaled by the trace of matrix "s" over the number of samples
              s / (sum(diag(s)) / nsamp)
            })

            ## eigen-decomposition
            eig <- eigen(genmat)

            eig$values <- eig$values[1:eigen.cnt]
            eig$vectors <- eig$vectors[,1:eigen.cnt]
            
            rownames(eig$vectors) <- seqGetData(gdsobj, "sample.id")
            
            names(eig) <- c("eigenval", "eigenvect")
            eig
          })
