
## applyMethod
setMethod("applyMethod",
          c(gdsobj="SeqVarGDSClass", FUN="function", variant="character"),
          function(gdsobj, FUN, variant, sample=NULL, ...) {
            .applyMethod(gdsobj, FUN, variant, sample, ...)
          })

setMethod("applyMethod",
          c(gdsobj="SeqVarGDSClass", FUN="function", variant="numeric"),
          function(gdsobj, FUN, variant, sample=NULL, ...) {
            .applyMethod(gdsobj, FUN, variant, sample, ...)
          })

setMethod("applyMethod",
          c(gdsobj="SeqVarGDSClass", FUN="function", variant="GRanges"),
          function(gdsobj, FUN, variant, sample=NULL, ...) {
            .applyMethod(gdsobj, FUN, .rangesToID(gdsobj, variant), sample, ...)
          })

setMethod("applyMethod",
          c(gdsobj="SeqVarGDSClass", FUN="function", variant="missing"),
          function(gdsobj, FUN, variant, sample=NULL, ...) {
            .applyMethod(gdsobj, FUN, variant.id=NULL, sample, ...)
          })


## allele methods
setMethod("refChar",
          "SeqVarGDSClass",
          function(gdsobj) {
            .parseRefAllele(seqGetData(gdsobj, "allele"))
          })

setMethod("altChar",
          "SeqVarGDSClass",
          function(gdsobj, n=0) {
            .parseAltAllele(seqGetData(gdsobj, "allele"), n=n)
          })

setMethod("nAlleles",
          "SeqVarGDSClass",
          function(gdsobj) {
            #.parseNumAlleles(seqGetData(gdsobj, "allele"))
            seqNumAllele(gdsobj)
          })


## descriptive methods
setMethod("isVariant",
          "SeqVarGDSClass",
          function(gdsobj, use.names=FALSE) {
            var <- seqApply(gdsobj, "genotype",
                            function(x) {colSums(x, na.rm=TRUE) > 0},
                            margin="by.variant", as.is="list")
            var <- matrix(unlist(var, use.names=FALSE), ncol=length(var),
                          dimnames=list(sample=NULL, variant=NULL))
            if (use.names) .applyNames(gdsobj, var) else var
          })
 
setMethod("isSNV",
          "SeqVarGDSClass",
          function(gdsobj, biallelic=TRUE) {
            a <- seqGetData(gdsobj, "allele")
            if (biallelic) {
              nchar(a) == 3
            } else {
              .maxAlleleLength(a) == 1
            }
          })


## data retrieval and formatting
setMethod("getGenotype",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            gc <- seqApply(gdsobj, c(geno="genotype", phase="phase"),
                           function(x) {sep=ifelse(x$phase, "|", "/")
                                        paste0(x$geno[1,], sep, x$geno[2,])},
                           margin="by.variant", as.is="list")
            gc <- matrix(unlist(gc, use.names=FALSE), ncol=length(gc),
                         dimnames=list(sample=NULL, variant=NULL))
            gc[gc == "NA/NA"] <- NA
            if (use.names) .applyNames(gdsobj, gc) else gc
          })

setMethod("getGenotypeAlleles",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE, sort=FALSE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            gc <- seqApply(gdsobj,
                           c(geno="genotype", phase="phase", allele="allele"),
                           function(x, sort) {
                             alleles <- unlist(strsplit(x$allele, ",", fixed=TRUE),
                                               use.names=FALSE)
                             names(alleles) <- 0:(length(alleles) - 1)
                             a <- alleles[as.character(x$geno[1,])]
                             b <- alleles[as.character(x$geno[2,])]
                             if (sort) {
                               paste(pmin(a,b), pmax(a,b), sep="/")
                             } else  {
                               sep=ifelse(x$phase, "|", "/")
                               paste0(a, sep, b)
                             }
                           }, margin="by.variant", as.is="list", sort=sort)
            gc <- matrix(unlist(gc, use.names=FALSE), ncol=length(gc),
                         dimnames=list(sample=NULL, variant=NULL))
            gc[gc == "NA/NA"] <- NA
            if (use.names) .applyNames(gdsobj, gc) else gc
          })

setMethod("refDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            rd <- seqApply(gdsobj, "genotype",
                           function(x) {colSums(x == 0)},
                           margin="by.variant", as.is="list")
            rd <- matrix(unlist(rd, use.names=FALSE), ncol=length(rd),
                         dimnames=list(sample=NULL, variant=NULL))
            if (use.names) .applyNames(gdsobj, rd) else rd
          })

setMethod("altDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            d <- seqApply(gdsobj, "genotype",
                           function(x) {colSums(x != 0)},
                           margin="by.variant", as.is="list")
            d <- matrix(unlist(d, use.names=FALSE), ncol=length(d),
                         dimnames=list(sample=NULL, variant=NULL))
            if (use.names) .applyNames(gdsobj, d) else d
          })

setMethod("alleleDosage",
          c("SeqVarGDSClass", "numeric"),
          function(gdsobj, n=0, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            if (length(n) == 1) n <- rep(n, .nVar(gdsobj))
            stopifnot(length(n) == .nVar(gdsobj))
            stopifnot(all(n <= nAlleles(gdsobj)))
            d <- seqApply(gdsobj, "genotype",
                          function(index, x) {colSums(x == n[index])},
                          margin="by.variant", as.is="list",
                          var.index="relative")
            d <- matrix(unlist(d, use.names=FALSE), ncol=length(d),
                         dimnames=list(sample=NULL, variant=NULL))
            if (use.names) .applyNames(gdsobj, d) else d
          })

## setMethod("alleleDosage",
##           c("SeqVarGDSClass", "list"),
##           function(gdsobj, n, use.names=TRUE) {
##             stopifnot(length(n) == .nVar(gdsobj))
##             d <- seqApply(gdsobj, "genotype",
##                           function(index, x) {
##                               tmp <- matrix(x %in% n[[index]], ncol=ncol(x),
##                                              nrow=nrow(x))
##                               tmp[is.na(x)] <- NA
##                               colSums(tmp)
##                           },
##                           margin="by.variant", as.is="list",
##                           var.index="relative")
##             d <- matrix(unlist(d, use.names=FALSE), ncol=length(d),
##                          dimnames=list(sample=NULL, variant=NULL))
##             if (use.names) .applyNames(gdsobj, d) else d
##           })

setMethod("alleleDosage",
          c("SeqVarGDSClass", "list"),
          function(gdsobj, n, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            stopifnot(length(n) == .nVar(gdsobj))
            samp.names <- if (use.names) seqGetData(gdsobj, "sample.id") else NULL
            d <- seqApply(gdsobj, "genotype",
                          function(index, x) {
                              tmp <- lapply(n[[index]], function(allele) colSums(x == allele))
                              matrix(unlist(tmp, use.names=FALSE), ncol=length(tmp),
                                     dimnames=list(sample=samp.names, allele=n[[index]]))
                          },
                          margin="by.variant", as.is="list",
                          var.index="relative")
            if (use.names) names(d) <- seqGetData(gdsobj, "variant.id")
            d
          })

setMethod("expandedAltDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            samp.names <- if (use.names) seqGetData(gdsobj, "sample.id") else NULL
            variant.id <- if (use.names) seqGetData(gdsobj, "variant.id") else NULL
            n <- nAlleles(gdsobj) - 1
            d <- seqApply(gdsobj, "genotype",
                          function(index, x) {
                              tmp <- lapply(1:n[index], function(allele) colSums(x == allele))
                              var.names <- rep(variant.id[index], n[index])
                              matrix(unlist(tmp, use.names=FALSE), ncol=length(tmp),
                                     dimnames=list(sample=samp.names, variant=var.names))
                          },
                          margin="by.variant", as.is="list",
                          var.index="relative")
            d <- do.call(cbind, d)
            if (use.names) names(dimnames(d)) <- c("sample", "variant")
            d
          })

setMethod("getVariableLengthData",
          c("SeqVarGDSClass", "character"),
          function(gdsobj, var.name, use.names=TRUE) {
            var.list <- seqApply(gdsobj, var.name, function(x) {x},
                                 margin="by.variant", as.is="list")
            .ncol <- function(x) ifelse(is.null(x), 0, ncol(x))
            n.alleles <- sapply(var.list, .ncol)
            nsamp <- .nSamp(gdsobj)
            if (length(unique(n.alleles)) > 1) {
                ## fill in missing values
                maxn <- max(n.alleles)
                var.list <- lapply(var.list, function(x) {
                    cbind(x, matrix(NA, nrow=nsamp, ncol=(maxn - .ncol(x))))
                })
            }
            var <- array(unlist(var.list, use.names=FALSE),
                         dim=c(nrow(var.list[[1]]), ncol(var.list[[1]]),
                             length(var.list)))
            var <- aperm(var, c(2,1,3))
            
            ## if first array dimension is 1, simplify to a matrix
            if (dim(var)[1] == 1) {
              var <- var[1,,]
            }
            if (length(dim(var)) == 3) {
              dimnames(var) <- list(n=NULL, sample=NULL, variant=NULL)
            } else {
              dimnames(var) <- list(sample=NULL, variant=NULL)
            }
            
            if (use.names) .applyNames(gdsobj, var) else var
          })

## metrics    
setMethod("titv",
          "SeqVarGDSClass",
          function(gdsobj, by.sample=FALSE, use.names=FALSE) {
            ref <- refChar(gdsobj)
            alt <- altChar(gdsobj)
            ti <- .isTransition(ref, alt)
            tv <- .isTransversion(ref, alt)
            if (by.sample) {
              isVar <- function(x) {colSums(x, na.rm=TRUE) > 0}
              tisum <- integer(.nSamp(gdsobj))
              tvsum <- integer(.nSamp(gdsobj))
              seqApply(gdsobj, "genotype",
                       function(index, x) {
                         tisum <<- tisum + (ti[index] & isVar(x));
                         tvsum <<- tvsum + (tv[index] & isVar(x))
                       },
                       margin="by.variant", as.is="none", var.index="relative")
              titv <- tisum / tvsum
              if (use.names)
                names(titv) <- seqGetData(gdsobj, "sample.id")
              titv
            } else {
              sum(ti) / sum(tv)
            }
          })

## n=0: REF allele freq
## n>0: freq of nth ALT allele
setMethod("alleleFrequency",
          "SeqVarGDSClass",
          function(gdsobj, n=0, use.names=FALSE) {
            ## af <- seqApply(gdsobj, "genotype",
            ##                function(x) {mean(x == n, na.rm=TRUE)},
            ##                margin="by.variant", as.is="double")
            af <- seqAlleleFreq(gdsobj, ref.allele=n)
            ## frequency of alt=2, etc. should be 0 if there is no such allele
            af[is.na(af)] <- 0
            if (use.names)
              names(af) <- seqGetData(gdsobj, "variant.id")
            af
          })

setMethod("missingGenotypeRate",
          "SeqVarGDSClass",
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE) {
            margin <- match.arg(margin)
            if (margin == "by.variant") {
              miss <- seqMissing(gdsobj, per.variant=TRUE)
              if (use.names)
                names(miss) <- seqGetData(gdsobj, "variant.id")
            } else {
              miss <- seqMissing(gdsobj, per.variant=FALSE)
              if (use.names)
                names(miss) <- seqGetData(gdsobj, "sample.id")
            } 
            miss

##             if (margin == "by.variant") {
##               miss <- seqApply(gdsobj, "genotype",
##                                function(x) {sum(is.na(x[1,]))},
##                                margin=margin, as.is="integer")
##               if (use.names)
##                 names(miss) <- seqGetData(gdsobj, "variant.id")
##               miss / .nSamp(gdsobj)
##             } else {
##               miss <- integer(.nSamp(gdsobj))
##               ## use "<<-" operator to find "miss" in the parent environment
##               seqApply(gdsobj, "genotype",
##                        function(x) {miss <<- miss + is.na(x[1,])},
##                        margin="by.variant", as.is="none")
##               if (use.names)
##                 names(miss) <- seqGetData(gdsobj, "sample.id")
##               miss / .nVar(gdsobj)
##             } 

## REPLACE when by.sample works in seqApply
##              miss <- seqApply(gdsobj, "genotype",
##                               function(x) {sum(is.na(x[1,]))},
##                               margin=margin, as.is="integer")
##             if (margin == "by.variant") {
##               miss / .nSamp(gdsobj)
##             } else {
##               miss / .nVar(gdsobj)
##             } 
          })

setMethod("heterozygosity",
          "SeqVarGDSClass",
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE) {
            margin <- match.arg(margin)
            if (margin == "by.variant") {
            het <- seqApply(gdsobj, "genotype",
                            function(x) {
                              sum(x[1,] != x[2,], na.rm=TRUE) /
                              sum(!is.na(x[1,]) & !is.na(x[2,]))
                            },
                            margin=margin, as.is="double")
            if (use.names)
              names(het) <- seqGetData(gdsobj, "variant.id")
            het
          } else {
            het <- integer(.nSamp(gdsobj))
            nonmiss <- integer(.nSamp(gdsobj))
            ## use "<<-" operator to find "het" in the parent environment
            seqApply(gdsobj, "genotype",
                     function(x) {
                       nm <- !is.na(x[1,]) & !is.na(x[2,])
                       het <<- het + (x[1,] != x[2,] & nm)
                       nonmiss <<- nonmiss + nm
                     },
                     margin="by.variant", as.is="none")
            if (use.names)
              names(het) <- seqGetData(gdsobj, "sample.id")
            het / nonmiss
          }
            
## REPLACE when by.sample works in seqApply
##             het <- seqApply(gdsobj, "genotype",
##                             function(x) {
##                               sum(x[1,] != x[2,], na.rm=TRUE) /
##                                 sum(!is.na(x[1,]) & !is.na(x[2,]))
##                             },
##                             margin=margin, as.is="integer")
          })

setMethod("homozygosity",
          "SeqVarGDSClass",
          function(gdsobj, allele=c("any", "ref", "alt"),
                   margin=c("by.variant", "by.sample"), use.names=FALSE) {
            hom.func <- switch(match.arg(allele),
                               any=function(a,b) {a == b},
                               ref=function(a,b) {a == b & a == 0},
                               alt=function(a,b) {a == b & a > 0})
            margin <- match.arg(margin)
            if (margin == "by.variant") {
              hom <- seqApply(gdsobj, "genotype",
                              function(x) {
                                sum(hom.func(x[1,], x[2,]), na.rm=TRUE) /
                                  sum(!is.na(x[1,]) & !is.na(x[2,]))
                              },
                              margin=margin, as.is="double")
              if (use.names)
                names(hom) <- seqGetData(gdsobj, "variant.id")
              hom
            } else {
              hom <- integer(.nSamp(gdsobj))
              nonmiss <- integer(.nSamp(gdsobj))
              ## use "<<-" operator to find "hom" in the parent environment
              seqApply(gdsobj, "genotype",
                       function(x) {
                         nm <- !is.na(x[1,]) & !is.na(x[2,])
                         hom <<- hom + (hom.func(x[1,], x[2,]) & nm)
                         nonmiss <<- nonmiss + nm
                       },
                       margin="by.variant", as.is="none")
              if (use.names)
                names(hom) <- seqGetData(gdsobj, "sample.id")
              hom / nonmiss
            }
          })


setMethod("hethom",
          "SeqVarGDSClass",
          function(gdsobj, use.names=FALSE) {
            hom.func <- function(a,b) {a == b & a > 0}
            het <- integer(.nSamp(gdsobj))
            hom <- integer(.nSamp(gdsobj))
            seqApply(gdsobj, "genotype",
                     function(x) {
                       nm <- !is.na(x[1,]) & !is.na(x[2,])
                       het <<- het + (x[1,] != x[2,] & nm)
                       hom <<- hom + (hom.func(x[1,], x[2,]) & nm)
                     },
                     margin="by.variant", as.is="none")
             if (use.names)
               names(het) <- seqGetData(gdsobj, "sample.id")
             het / hom
           })
            

setMethod("meanBySample",
          "SeqVarGDSClass",
          function(gdsobj, var.name, use.names=FALSE) {
            ns <- .nSamp(gdsobj)
            tot <- double(ns)
            nm <- double(ns)
            seqApply(gdsobj, var.name,
                     function(x) {
                       val <- !is.na(x)
                       tot[val] <<- tot[val] + x[val]
                       nm[val] <<- nm[val] + 1
                     }, margin="by.variant", as.is="none")
            if (use.names)
              names(tot) <- seqGetData(gdsobj, "sample.id")
            tot / nm
          })


setMethod("countSingletons",
          "SeqVarGDSClass",
          function(gdsobj, use.names=FALSE) {
            st <- integer(.nSamp(gdsobj))
            seqApply(gdsobj, "genotype",
                     function(x) {
                       nonref <- which(colSums(x, na.rm=TRUE) > 0)
                       if (length(nonref) == 1) {
                           st[nonref] <<- st[nonref] + 1
                       }
                     },
                     margin="by.variant", as.is="none")
             if (use.names)
               names(st) <- seqGetData(gdsobj, "sample.id")
             st
           })
      
