
setMethod("refFrac",
          "SeqVarGDSClass",
          function(gdsobj, use.names=TRUE) {
            vars <- seqSummary(gdsobj, "annotation/format", check="none", verbose=FALSE)$ID
            if (!("AD" %in% vars)) {
              stop("annotation/format/AD must be present to compute allelic balance")
            }
            
            rf <- seqApply(gdsobj, "annotation/format/AD",
                           function(x) {x[,1] / rowSums(x)},
                           margin="by.variant", as.is="list")
            rf <- matrix(unlist(rf), ncol=length(rf),
                         dimnames=list(sample=NULL, variant=NULL))
            if (use.names) .applyNames(gdsobj, rf) else rf
          })


setMethod("refFracOverHets",
          "SeqVarGDSClass",
          function(gdsobj, FUN=mean, use.names=TRUE) {
            vars <- seqSummary(gdsobj, "annotation/format", check="none", verbose=FALSE)$ID
            if (!("AD" %in% vars)) {
              stop("annotation/format/AD must be present to compute allelic balance")
            }

            rf <- seqApply(gdsobj,
                           c(ad="annotation/format/AD", geno="genotype"),
                           function(x) {
                             het <- x$geno[1,] != x$geno[2,]
                             FUN(x$ad[het,1,drop=FALSE] /
                                 rowSums(x$ad[het,,drop=FALSE]), na.rm=TRUE)
                           },
                           margin="by.variant", as.is="double")   
            if (use.names)
              names(rf) <- seqGetData(gdsobj, "variant.id")
            rf
          })


## code by genotype
.callType <- function(x) {
  pcol <- character(length(x))
  a <- substr(x, 1, 1)
  b <- substr(x, 3, 3)
  pcol[a == "0" | b == "0"] <- "refhet"
  pcol[a == "0" & b == "0"] <- "refhom"
  pcol[a != "0" & b != "0" & a == b] <- "althom"
  pcol[a != "0" & b != "0" & a != b] <- "althet"
  pcol[is.na(x)] <- "na"
  pcol
}


## title with variant.id, chr:position, SNV status, genotype counts
.varTitle <- function(gdsobj) {
  var.id <- seqGetData(gdsobj, "variant.id")
  loc <- paste(seqGetData(gdsobj, "chromosome"),
               seqGetData(gdsobj, "position"), sep=":")
  snv <- isSNV(gdsobj, biallelic=FALSE)
  nalleles <- nAlleles(gdsobj)
  counts <- .countGenotypes(gdsobj)
  paste0("Variant ", var.id, ", chr", loc, ", ",
        ifelse(snv, "SNV", "INDEL"), ", ",
        "nAlt=", nalleles-1, "\n",
         apply(counts, 1, function(x) paste(rev(paste(names(counts), x, sep="=")), collapse=", ")))
}


setMethod("refFracPlot",
          "SeqVarGDSClass",
          function(gdsobj, variant.id, highlight=NULL,...) {
            filt.orig <- seqGetFilter(gdsobj)$variant.sel
            seqSetFilter(gdsobj, variant.id=variant.id, verbose=FALSE)

            ## number of reference reads
            num.ref <- seqApply(gdsobj, "annotation/format/AD",
                           function(x) {x[,1]},
                           margin="by.variant", as.is="list")
            num.ref <- matrix(unlist(num.ref), ncol=length(num.ref))

            ## total read depth
            tot.reads <- seqApply(gdsobj, "annotation/format/AD",
                           rowSums,
                           margin="by.variant", as.is="list")
            tot.reads <- matrix(unlist(tot.reads), ncol=length(tot.reads))

            ## reference allele fraction
            ref.frac <- num.ref / tot.reads

            ## color-code by genotype call
            geno <- getGenotype(gdsobj)
            dimnames(ref.frac) <- dimnames(geno)
            dimnames(tot.reads) <- dimnames(geno)

            ## title with variant.id, chr:position, SNV status
            title <- .varTitle(gdsobj)

            ## plot in same order as variant.id
            for (i in match(variant.id, seqGetData(gdsobj, "variant.id"))) {            
              ## median fraction for hets
              calls <- .callType(geno[,i])
              refhet <- calls == "refhet"
              med.frac <- median(ref.frac[refhet, i], na.rm=TRUE)
              ## hets significantly different from 0.5
              pch <- rep(16, length(refhet))
              alpha <- 0.05 / sum(refhet)  # Bonferroni correction
              for (j in which(refhet)) {
                if (tot.reads[j,i] > 0) {
                  sig <- binom.test(num.ref[j,i], tot.reads[j,i], 0.5,
                                    alternative="two.sided")$p.value
                  if (sig < alpha) pch[j] <- 17
                }
              }
              ## missing genotype
              pch[calls == "na"] <- 4

              colmap <- c("refhom"="#1b9e77",
                          "refhet"="#d95f02",
                          "althom"="#7570b3",
                          "althet"="#e7298a",
                          "na"="black")
              
              plot(ref.frac[,i], tot.reads[,i],
                   col=adjustcolor(colmap[calls], alpha.f=0.5), pch=pch,
                   xlab="Fraction of reference allele reads",
                   ylab="Total reads (unfiltered)",
                   main=title[i], xlim=c(0,1), ...)
              abline(v=0.5)
              abline(v=med.frac, col=colmap["refhet"])
              
              ## highlight samples?
              if (!is.null(highlight)) {
                  text(ref.frac[highlight[[i]],i], tot.reads[highlight[[i]],i],
                       col="black")
              }
            }
            
            seqSetFilter(gdsobj, variant.sel=filt.orig, verbose=FALSE)
          })

