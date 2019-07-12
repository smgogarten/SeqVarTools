## constructor
SeqVarData <- function(gds, sampleData, variantData) {
    closeOnFail <- FALSE
    if (is.character(gds) & length(gds) == 1) {
        gds <- seqOpen(gds)
        closeOnFail <- TRUE
    } else if (!is(gds, "SeqVarGDSClass")) {
        stop("gds must be a SeqVarGDSClass object or a valid filename")
    }
    class(gds) <- "SeqVarData"

    if (missing(sampleData)) {
        sampleData <- AnnotatedDataFrame(data.frame(matrix(nrow=.nSampUnfiltered(gds), ncol=0)))
    }
    gds@sampleData <- sampleData
    
    if (missing(variantData)) {
        variantData <- AnnotatedDataFrame(data.frame(matrix(nrow=.nVarUnfiltered(gds), ncol=0)))
    }
    gds@variantData <- variantData
    
    check <- validObject(gds, test=TRUE)
    if (is.character(check)) {
        if (closeOnFail) seqClose(gds)
        stop(check)
    }
    gds
}

setValidity("SeqVarData",
            function(object) {
                if (ncol(sampleData(object)) > 0) {
                    if (!("sample.id" %in% varLabels(sampleData(object)))) {
                        return("sampleData does not have column sample.id")
                    }
                    if (!identical(seqGetData(object, "sample.id"),
                                   sampleData(object)$sample.id)) {
                        return("sample.id in sampleData does not match GDS file")
                    }
                }
                if (ncol(variantData(object)) > 0) {
                    if (!("variant.id" %in% varLabels(variantData(object)))) {
                        return("variantData does not have column variant.id")
                    }
                    if (!identical(seqGetData(object, "variant.id"),
                                   variantData(object)$variant.id)) {
                        return("variant.id in variantData does not match GDS file")
                    }
                }
                TRUE
            })

setMethod("show",
          "SeqVarData",
          function(object) {
              cat(class(object), "object\n")
              cat(" | GDS:\n")
              print(object)
              cat(" | sampleData:\n")
              show(sampleData(object))
              cat(" | variantData:\n")
              show(variantData(object))
          })


setMethod("sampleData",
          "SeqVarData",
          function(x) {
              sd <- x@sampleData
              samp.sel <- seqGetFilter(x)$sample.sel
              if (sum(samp.sel) == length(samp.sel)) {
                  return(sd)
              } else {
                  return(sd[samp.sel,])
              }
          })

setReplaceMethod("sampleData",
                 c("SeqVarData", "AnnotatedDataFrame"),
                 function(x, value) {
                     seqSetFilter(x, sample.sel=1:.nSampUnfiltered(x), action="push+set", verbose=FALSE)
                     if (!identical(seqGetData(x, "sample.id"), value$sample.id)) {
                         stop("sample.id does not match GDS file")
                     }
                     seqSetFilter(x, action="pop", verbose=FALSE)
                     slot(x, "sampleData") <- value
                     x
                 })

            
setMethod("variantData",
          "SeqVarData",
          function(x) {
              vd <- x@variantData
              var.sel <- seqGetFilter(x)$variant.sel
              if (sum(var.sel) == length(var.sel)) {
                  return(vd)
              } else {
                  return(vd[var.sel,])
              }
          })

setReplaceMethod("variantData",
                 c("SeqVarData", "AnnotatedDataFrame"),
                 function(x, value) {
                     seqSetFilter(x, variant.sel=1:.nVarUnfiltered(x), action="push+set", verbose=FALSE)
                     if (!identical(seqGetData(x, "variant.id"), value$variant.id)) {
                         stop("variant.id does not match GDS file")
                     }
                     seqSetFilter(x, action="pop", verbose=FALSE)
                     slot(x, "variantData") <- value
                     x
                 })


setMethod("granges",
          "SeqVarData",
          function(x, ...) {
              gr <- callNextMethod(x, ...)
              df <- pData(variantData(x))
              df$variant.id <- NULL # redundant with names of gr
              mcols(gr) <- cbind(mcols(gr), as(df, "DataFrame"))
              gr
          })


setMethod("validateSex",
          "SeqVarData",
          function(x) {
              sex <- sampleData(x)$sex
              if (!is.null(sex)) {          
                  if (all(sex %in% c(1,2,NA))) {
                      sex <- c("M", "F")[sex]
                  }
                  if (!all(sex %in% c("M", "F", NA))) {
                      sex <- NULL
                  }
              }
              sex
          })


setMethod("alleleFrequency",
          "SeqVarData",
          function(gdsobj, n=0, use.names=FALSE, sex.adjust=TRUE, genome.build=c("hg19", "hg38")) {
              if (!sex.adjust) {
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }
              
              # check chromosome
              chr <- chromWithPAR(gdsobj, genome.build)
              if (!any(chr %in% c("X", "Y"))) {
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }

              # check sex
              sex <- validateSex(gdsobj)
              if (is.null(sex)) {
                  warning("No valid sex coding provided in sampleData. Frequencies will not be calculated correctly for X and Y chromosomes.")
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }
              female <- sex %in% "F"
              male <- sex %in% "M"

              # empty vector to fill in
              freq <- rep(NA, length(chr))
              if (use.names) names(freq) <- seqGetData(gdsobj, "variant.id")

              # autosomal
              auto <- !(chr %in% c("X", "Y"))
              if (any(auto)) {
                  seqSetFilter(gdsobj, variant.sel=auto, action="push+intersect", verbose=FALSE)
                  freq[auto] <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }
              
              # X chrom
              X <- chr %in% "X"
              if (any(X)) {
                  seqSetFilter(gdsobj, sample.sel=female, variant.sel=X, action="push+intersect", verbose=FALSE)
                  freq.X.F <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  n.F <- .nSampObserved(gdsobj)
                  count.X.F <- freq.X.F * 2*n.F
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=X, action="push+intersect", verbose=FALSE)
                  freq.X.M <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  n.M <- .nSampObserved(gdsobj)
                  count.X.M <- freq.X.M * n.M
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
                  
                  freq[X] <- (count.X.F + count.X.M)/(2*n.F + n.M)
              }

              # Y chrom
              Y <- chr %in% "Y"
              if (any(Y)) {
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=Y, action="push+intersect", verbose=FALSE)
                  freq[Y] <- callNextMethod(gdsobj, n=n, use.names=use.names)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }

              freq
          })


setMethod("alleleCount",
          "SeqVarData",
          function(gdsobj, n=0, use.names=FALSE, sex.adjust=TRUE, genome.build=c("hg19", "hg38")) {
              if (!sex.adjust) {
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }
              
              # check chromosome
              chr <- chromWithPAR(gdsobj, genome.build)
              if (!any(chr %in% c("X", "Y"))) {
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }

              # check sex
              sex <- validateSex(gdsobj)
              if (is.null(sex)) {
                  warning("No valid sex coding provided in sampleData. Counts will not be calculated correctly for X and Y chromosomes.")
                  return(callNextMethod(gdsobj, n=n, use.names=use.names))
              }
              female <- sex %in% "F"
              male <- sex %in% "M"

              # empty vector to fill in
              count <- rep(NA, length(chr))
              if (use.names) names(count) <- seqGetData(gdsobj, "variant.id")

              # autosomal
              auto <- !(chr %in% c("X", "Y"))
              if (any(auto)) {
                  seqSetFilter(gdsobj, variant.sel=auto, action="push+intersect", verbose=FALSE)
                  count[auto] <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }
              
              # X chrom
              X <- chr %in% "X"
              if (any(X)) {
                  seqSetFilter(gdsobj, sample.sel=female, variant.sel=X, action="push+intersect", verbose=FALSE)
                  count.X.F <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=X, action="push+intersect", verbose=FALSE)
                  count.X.M <- callNextMethod(gdsobj, n=n, use.names=FALSE)
                  count.X.M <- count.X.M / 2 ## correct for X chromosome coded as diploid in males
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
                  
                  count[X] <- count.X.F + count.X.M
              }

              # Y chrom
              Y <- chr %in% "Y"
              if (any(Y)) {
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=Y, action="push+intersect", verbose=FALSE)
                  count[Y] <- callNextMethod(gdsobj, n=n, use.names=use.names)
                  count[Y] <- count[Y] / 2 ## correct for Y chromosome coded as diploid in males
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }

              count
          })


setMethod("minorAlleleCount",
          "SeqVarData",
          function(gdsobj, use.names=FALSE, sex.adjust=TRUE, genome.build=c("hg19", "hg38")) {
              mac <- function(gdsobj, use.names=use.names) {
                  ref.cnt <- alleleCount(gdsobj, n=0, use.names=use.names)
                  n.obs <- .nSampObserved(gdsobj)
                  pmin(ref.cnt, 2*n.obs - ref.cnt)
              }
              
              if (!sex.adjust) {
                  return(mac(gdsobj, use.names=use.names))
              }
              
              # check chromosome
              chr <- chromWithPAR(gdsobj, genome.build)
              if (!any(chr %in% c("X", "Y"))) {
                  return(mac(gdsobj, use.names=use.names))
              }

              # check sex
              sex <- validateSex(gdsobj)
              if (is.null(sex)) {
                  warning("No valid sex coding provided in sampleData. Counts will not be calculated correctly for X and Y chromosomes.")
                  return(mac(gdsobj, use.names=use.names))
              }
              female <- sex %in% "F"
              male <- sex %in% "M"

              # empty vector to fill in
              count <- rep(NA, length(chr))
              possible <- rep(NA, length(chr))
              if (use.names) names(count) <- seqGetData(gdsobj, "variant.id")

              # autosomal
              auto <- !(chr %in% c("X", "Y"))
              if (any(auto)) {
                  seqSetFilter(gdsobj, variant.sel=auto, action="push+intersect", verbose=FALSE)
                  count[auto] <- seqAlleleCount(gdsobj)
                  possible[auto] <- 2 * .nSampObserved(gdsobj)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }
              
              # X chrom
              X <- chr %in% "X"
              if (any(X)) {
                  seqSetFilter(gdsobj, sample.sel=female, variant.sel=X, action="push+intersect", verbose=FALSE)
                  count.X.F <- seqAlleleCount(gdsobj)
                  possible.X.F <- 2 * .nSampObserved(gdsobj)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=X, action="push+intersect", verbose=FALSE)
                  count.X.M <- seqAlleleCount(gdsobj)
                  count.X.M <- count.X.M / 2 ## correct for X chromosome coded as diploid in males
                  possible.X.M <- .nSampObserved(gdsobj)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
                  
                  count[X] <- count.X.F + count.X.M
                  possible[X] <- possible.X.F + possible.X.M
              }

              # Y chrom
              Y <- chr %in% "Y"
              if (any(Y)) {
                  seqSetFilter(gdsobj, sample.sel=male, variant.sel=Y, action="push+intersect", verbose=FALSE)
                  count[Y] <- seqAlleleCount(gdsobj)
                  count[Y] <- count[Y] / 2 ## correct for Y chromosome coded as diploid in males
                  possible[Y] <- .nSampObserved(gdsobj)
                  seqSetFilter(gdsobj, action="pop", verbose=FALSE)
              }

              pmin(count, possible - count)
          })
