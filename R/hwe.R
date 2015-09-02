
setMethod("hwe",
          "SeqVarGDSClass",
          function(gdsobj, use.names=FALSE) {
            ## genotype counts
            nAA <- seqApply(gdsobj, "genotype",
                            function(x) {sum(x[1,] == 0 & x[2,] == 0, na.rm=TRUE)},
                            margin="by.variant", as.is="integer")
            nAa <- seqApply(gdsobj, "genotype",
                            function(x) {sum((x[1,] == 0 & x[2,] == 1) |
                                             (x[1,] == 1 & x[2,] == 0), na.rm=TRUE)},
                            margin="by.variant", as.is="integer")
            naa <- seqApply(gdsobj, "genotype",
                            function(x) {sum(x[1,] == 1 & x[2,] == 1, na.rm=TRUE)},
                            margin="by.variant", as.is="integer")
            pv <- HWExact(data.frame(nAA, nAa, naa))

            ## set non-biallelic variants to NA
            pv[nAlleles(gdsobj) != 2] <- NA

            ## set pvalue to NA for monomorphic SNPs
            pv[(nAa + naa) == 0 | (nAa + nAA) == 0] <- NA

            if (use.names)
              names(pv) <- seqGetData(gdsobj, "variant.id")
            pv
          })
          
