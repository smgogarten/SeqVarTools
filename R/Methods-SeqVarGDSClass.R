
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
          function(gdsobj, use.names=FALSE, parallel=FALSE) {
            var <- seqApply(gdsobj, "genotype",
                            function(x) {colSums(x, na.rm=TRUE) > 0},
                            margin="by.variant", as.is="list",
                            parallel=parallel)
            var <- matrix(unlist(var, use.names=FALSE), ncol=length(var))
            if (use.names) {
                dimnames(var) <- list(sample=NULL, variant=NULL)
                .applyNames(gdsobj, var)
            }else var
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
          function(gdsobj, use.names=TRUE, parallel=FALSE) {
              if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
              if (.ploidy(gdsobj) == 1) {
                  gc <- seqGetData(gdsobj, "genotype")[1,,]
              } else {
                  gc <- seqApply(gdsobj, c(geno="genotype", phase="phase"),
                                 function(x) {sep=ifelse(x$phase, "|", "/")
                                     paste0(x$geno[1,], sep, x$geno[2,])},
                                 margin="by.variant", as.is="list",
                                 parallel=parallel)
                  gc <- matrix(unlist(gc, use.names=FALSE), ncol=length(gc))
                  ## gc <- seqBlockApply(gdsobj, c(geno="genotype", phase="phase"),
                  ##                     function(x) {
                  ##                         sep <- ifelse(x$phase, "|", "/")
                  ##                         matrix(paste0(x$geno[1,,], sep, x$geno[2,,]), nrow=ncol(x$geno))
                  ##                     },
                  ##                margin="by.variant", as.is="list", ...)
                  ## gc <- do.call(cbind, gc)
                  gc[gc == "NA/NA"] <- NA
              }
              if (use.names) {
                  dimnames(gc) <- list(sample=NULL, variant=NULL)
                  .applyNames(gdsobj, gc)
              } else gc
          })

setMethod("getGenotypeAlleles",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE, sort=FALSE, parallel=FALSE) {
              if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
              if (.ploidy(gdsobj) == 1) {
                  gc <- seqApply(gdsobj,
                                 c(geno="genotype", allele="allele"),
                                 function(x) {
                                     alleles <- unlist(strsplit(x$allele, ",", fixed=TRUE),
                                                       use.names=FALSE)
                                     names(alleles) <- 0:(length(alleles) - 1)
                                     alleles[as.character(x$geno[1,])]
                                 }, margin="by.variant", as.is="list",
                                 parallel=parallel)
                  gc <- matrix(unlist(gc, use.names=FALSE), ncol=length(gc))
              } else {
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
                                 }, margin="by.variant", as.is="list", sort=sort,
                                 parallel=parallel)
                  gc <- matrix(unlist(gc, use.names=FALSE), ncol=length(gc))
                  ## gc <- seqBlockApply(gdsobj,
                  ##                c(geno="genotype", phase="phase", allele="allele"),
                  ##                function(x, sort) {
                  ##                    alleles <- strsplit(x$allele, ",", fixed=TRUE)
                  ##                    allele.map <- function(g, allele) {
                  ##                        names(allele) <- 0:(length(allele) - 1)
                  ##                        unname(allele[as.character(g)])
                  ##                    }
                  ##                    a <- mapply(allele.map, as.data.frame(x$geno[1,,]), alleles, USE.NAMES=FALSE)
                  ##                    b <- mapply(allele.map, as.data.frame(x$geno[2,,]), alleles, USE.NAMES=FALSE)
                  
                  ##                    if (sort) {
                  ##                        g <- paste(pmin(a,b), pmax(a,b), sep="/")
                  ##                    } else  {
                  ##                        sep=ifelse(x$phase, "|", "/")
                  ##                        g <- paste0(a, sep, b)
                  ##                    }
                  ##                    matrix(g, nrow=ncol(x$geno))
                  ##                }, margin="by.variant", as.is="list", sort=sort, ...)
                  ## gc <- do.call(cbind, gc)
                  gc[gc == "NA/NA"] <- NA
              }
              if (use.names) {
                  dimnames(gc) <- list(sample=NULL, variant=NULL)
                  .applyNames(gdsobj, gc)
              } else gc
          })

setMethod("refDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE, ...) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            ## d <- seqBlockApply(gdsobj, "genotype",
            ##                    function(x) {colSums(x == 0)},
            ##                    margin="by.variant", as.is="list", ...)
            ## d <- do.call(cbind, d)
            d <- seqGetData(gdsobj, "$dosage", ...)
            if (use.names) {
                #dimnames(d) <- list(sample=NULL, variant=NULL)
                d <- .applyNames(gdsobj, d)
            } else {
                dimnames(d) <- NULL
            }
            d
          })

setMethod("altDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE, sparse=FALSE, parallel=FALSE, ...) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            d <- seqBlockApply(gdsobj, "genotype",
                               function(x) {
                                   m <- colSums(x != 0)
                                   if (sparse) m <- Matrix(m, sparse=TRUE)
                                   m
                               }, margin="by.variant", as.is="list",
                               parallel=parallel, ...)
            #d <- do.call(cbind, d)
            d <- Reduce(cbind, d[-1], d[[1]])
            if (use.names) {
                dimnames(d) <- list(sample=NULL, variant=NULL)
                .applyNames(gdsobj, d)
            } else d
          })

## uses too much memory
## setMethod("altDosage2",
##           "SeqVarGDSClass",
##           function(gdsobj, use.names=TRUE, sparse=FALSE, ...) {
##             if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
##             d <- seqGetData(gdsobj, "$dosage_alt", ...)
##             if (sparse) d <- Matrix(d, sparse=TRUE)
##             if (use.names) {
##                 d <- .applyNames(gdsobj, d)
##             } else {
##                 dimnames(d) <- if (sparse) list(NULL,NULL) else NULL
##             }
##             d
##           })

setMethod("alleleDosage",
          c("SeqVarGDSClass", "numeric"),
          function(gdsobj, n=0, use.names=TRUE, parallel=FALSE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))
            if (length(n) == 1) n <- rep(n, .nVar(gdsobj))
            stopifnot(length(n) == .nVar(gdsobj))
            stopifnot(all(n <= nAlleles(gdsobj)))
            d <- seqApply(gdsobj, "genotype",
                          function(index, x) {colSums(x == n[index])},
                          margin="by.variant", as.is="list",
                          var.index="relative", parallel=parallel)
            d <- matrix(unlist(d, use.names=FALSE), ncol=length(d))
            if (use.names) {
                dimnames(d) <- list(sample=NULL, variant=NULL)
                .applyNames(gdsobj, d)
            } else d
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
          function(gdsobj, n, use.names=TRUE, parallel=FALSE) {
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
                          var.index="relative", parallel=parallel)
            if (use.names) names(d) <- seqGetData(gdsobj, "variant.id")
            d
          })

setMethod("expandedAltDosage",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE, sparse=FALSE, parallel=FALSE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))

            n <- nAlleles(gdsobj) - 1
            # if we have only biallelic variants, faster to use altDosage
            if (all(n == 1)) return(altDosage(gdsobj, use.names=use.names, sparse=sparse))
            
            samp.names <- if (use.names) seqGetData(gdsobj, "sample.id") else NULL
            variant.id <- if (use.names) seqGetData(gdsobj, "variant.id") else NULL
            d <- seqApply(gdsobj, "genotype",
                          function(index, x) {
                              tmp <- lapply(1:n[index], function(allele) colSums(x == allele))
                              var.names <- rep(variant.id[index], n[index])
                              m <- matrix(unlist(tmp, use.names=FALSE), ncol=length(tmp),
                                          dimnames=list(sample=samp.names, variant=var.names))
                              if (sparse) m <- Matrix(m, sparse=TRUE)
                              m
                          },
                          margin="by.variant", as.is="list",
                          var.index="relative", parallel=parallel)
            #https://stackoverflow.com/questions/37581417/getting-node-stack-overflow-when-cbind-multiple-sparse-matrices
            #d <- do.call(cbind, d)
            d <- Reduce(cbind, d[-1], d[[1]])
            if (use.names) names(dimnames(d)) <- c("sample", "variant")
            d
          })

setMethod("expandedVariantIndex",
          "SeqVarGDSClass",
          function(gdsobj) {
              nv <- .nVar(gdsobj)
              if (nv == 0) return(integer())
              nAlt <- nAlleles(gdsobj) - 1
              # keep variants with no alternate alleles
              nAlt[nAlt == 0] <- 1
              rep(1:nv, nAlt)
          })


### tidyverse version
## setMethod("variantInfo",
##           "SeqVarGDSClass",
##           function(gdsobj, alleles=TRUE, expanded=FALSE) {
##               x <- data.frame(variant.id=seqGetData(gdsobj, "variant.id"),
##                               chr=seqGetData(gdsobj, "chromosome"),
##                               pos=seqGetData(gdsobj, "position"),
##                               stringsAsFactors=FALSE)
##               if (nrow(x) == 0) return(x)
##               if (alleles) {
##                   allele <- seqGetData(gdsobj, "allele")
##                   x$ref <- .parseRefAllele(allele)
##                   x$alt <- .parseAltAllele(allele)
##               }
##               if (expanded) {
##                   if (alleles) {
##                       x <- separate_rows(x, .data$alt, sep=",")
##                   } else {
##                       x <- x[expandedVariantIndex(gdsobj),]
##                   }
##                   x <- group_by(x, .data$variant.id)
##                   x <- mutate(x, allele.index = 1:n())
##                   x <- as.data.frame(x)
##               }
##               x
##           })

setMethod("variantInfo",
          "SeqVarGDSClass",
          function(gdsobj, alleles=TRUE, expanded=FALSE) {
              x <- data.frame(variant.id=seqGetData(gdsobj, "variant.id"),
                              chr=seqGetData(gdsobj, "chromosome"),
                              pos=seqGetData(gdsobj, "position"),
                              stringsAsFactors=FALSE)
              if (nrow(x) == 0) return(x)
              if (alleles) {
                  allele <- seqGetData(gdsobj, "allele")
                  x <- as.data.table(x)
                  x[, `:=`(ref = .parseRefAllele(allele),
                           alt = .parseAltAllele(allele)
                           )]
              }
              if (expanded) {
                  x <- as.data.table(x)
                  if (alleles) {
                      alt <- NULL
                      x <- x[, strsplit(alt, ",", fixed=TRUE), by = c("variant.id", "chr", "pos", "ref", "alt")
                        ][, alt := NULL
                          ][, setnames(.SD, "V1", "alt")]
                  } else {
                      x <- x[expandedVariantIndex(gdsobj),]
                  }
                  x[, "allele.index" := seq_len(.N), by="variant.id"]
              }
              as.data.frame(x)
          })

setMethod("getVariableLengthData",
          c("SeqVarGDSClass", "character"),
          function(gdsobj, var.name, use.names=TRUE, parallel=FALSE) {
            var.list <- seqApply(gdsobj, var.name, function(x) {x},
                                 margin="by.variant", as.is="list",
                                 parallel=parallel)
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


setMethod("imputedDosage",
          "SeqVarGDSClass",
          function(gdsobj, dosage.field="DS", use.names=TRUE) {
            if (.emptyDim(gdsobj)) return(.emptyGenoMatrix(gdsobj, use.names=use.names))

            d <- seqGetData(gdsobj, paste0("annotation/format/", dosage.field))
            if (!is.null(names(d))) {
                if (!all(d$length == 1)) stop("multiple dosage values per variant")
                d <- d$data
            }
            if (use.names) {
                dimnames(d) <- list(sample=NULL, variant=NULL)
                .applyNames(gdsobj, d)
            } else d
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
          function(gdsobj, n=0, use.names=FALSE, parallel=FALSE) {
            ## af <- seqApply(gdsobj, "genotype",
            ##                function(x) {mean(x == n, na.rm=TRUE)},
            ##                margin="by.variant", as.is="double")
            af <- seqAlleleFreq(gdsobj, ref.allele=n, parallel=parallel)
            ## frequency of alt=2, etc. should be 0 if there is no such allele
            af[is.na(af)] <- 0
            if (use.names)
              names(af) <- seqGetData(gdsobj, "variant.id")
            af
          })

## n=0: REF allele count
## n>0: count of nth ALT allele
setMethod("alleleCount",
          "SeqVarGDSClass",
          function(gdsobj, n=0, use.names=FALSE, parallel=FALSE) {
            ac <- seqAlleleCount(gdsobj, ref.allele=n, parallel=parallel)
            ## count of alt=2, etc. should be 0 if there is no such allele
            ac[is.na(ac)] <- 0L
            if (use.names)
              names(ac) <- seqGetData(gdsobj, "variant.id")
            ac
          })

setMethod("minorAlleleCount",
          "SeqVarGDSClass",
          function(gdsobj, use.names=FALSE, parallel=FALSE) {
              ref.cnt <- seqAlleleCount(gdsobj, parallel=parallel)
              n.obs <- .nSampObserved(gdsobj)
              ac <- pmin(ref.cnt, 2L*n.obs - ref.cnt)
              if (use.names)
                  names(ac) <- seqGetData(gdsobj, "variant.id")
              ac
          })

setMethod("missingGenotypeRate",
          "SeqVarGDSClass",
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE, parallel=FALSE) {
            margin <- match.arg(margin)
            if (margin == "by.variant") {
              miss <- seqMissing(gdsobj, per.variant=TRUE, parallel=parallel)
              if (use.names)
                names(miss) <- seqGetData(gdsobj, "variant.id")
            } else {
              miss <- seqMissing(gdsobj, per.variant=FALSE, parallel=parallel)
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
          function(gdsobj, margin=c("by.variant", "by.sample"), use.names=FALSE, parallel=FALSE) {
            margin <- match.arg(margin)
            if (margin == "by.variant") {
            het <- seqApply(gdsobj, "genotype",
                            function(x) {
                              sum(x[1,] != x[2,], na.rm=TRUE) /
                              sum(!is.na(x[1,]) & !is.na(x[2,]))
                            },
                            margin=margin, as.is="double", parallel=parallel)
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
                   margin=c("by.variant", "by.sample"), use.names=FALSE, parallel=FALSE) {
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
                              margin=margin, as.is="double", parallel=parallel)
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
