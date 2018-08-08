.samplePairs1 <- function(samples) {
  stopifnot(all(c("sample.id", "subject.id") %in% names(samples)))

  ## keep only subjects with duplicate samples
  subjects <- unique(samples$subject.id[duplicated(samples$subject.id)])
  if (length(subjects) == 0) {
    stop("No duplicate samples - all values of subject.id are unique")
  }
  samples <- samples[samples$subject.id %in% subjects,]

  sample1 <- character(length(subjects))
  sample2 <- character(length(subjects))
  for (s in 1:length(subjects)) {
    samp.s <- samples$sample.id[samples$subject.id %in% subjects[s]]
    if (length(samp.s) > 2) {
      warning("More than two samples for subject ", s,
              "\nSelecting first two samples: ", samp.s[1], ", ", samp.s[2])
      samp.s <- samp.s[1:2]
    }
    sample1[s] <- samp.s[1]
    sample2[s] <- samp.s[2]
  }
  data.frame(sample1, sample2, row.names=subjects, stringsAsFactors=FALSE)
}

.samplePairs <- function(samples, all.pairs=TRUE) {
  stopifnot(all(c("sample.id", "subject.id") %in% names(samples)))

  ## keep only subjects with duplicate samples
  subjects <- unique(na.omit(samples$subject.id[duplicated(samples$subject.id)]))
  if (length(subjects) == 0) {
    stop("No duplicate samples - all values of subject.id are unique")
  }
  samples <- samples[samples$subject.id %in% subjects,]
  samples$order <- 1:nrow(samples)

  x <- merge(samples, samples, by="subject.id", suffixes = c(".1",".2"))
  x <- x[x$order.1 < x$order.2, c("subject.id", "sample.id.1", "sample.id.2")] 

  if (!all.pairs) {
      x <- x[!duplicated(x$subject.id),]
  }
  row.names(x) <- 1:nrow(x)
  return(x)
}

.genoMatch <- function(geno, check.phase) {
  allele1Match <- geno[1,1,] == geno[1,2,]
  allele2Match <- geno[2,1,] == geno[2,2,]
  if (check.phase) {
    allele1Match & allele2Match
  } else {
    allele12Match <- geno[1,1,] == geno[2,2,]
    allele21Match <- geno[2,1,] == geno[1,2,]
    (allele1Match & allele2Match) |
    (allele12Match & allele21Match)
  }
}

## setMethod("duplicateDiscordance",
##           c("SeqVarData", "missing"),
duplicateDiscordance1 <-
        function(gdsobj, match.samples.on="subject.id", check.phase=FALSE, verbose=TRUE) {
            ## samples should have columns of sample.id, subject.id
            ## find matching sample pairs for subjects (one sample per subject)
            samples <- pData(sampleData(gdsobj))
            if (!(match.samples.on %in% names(samples))) stop(sprintf("%s is not a column in sampleData", match.samples.on))
            samples <- samples[, c("sample.id", match.samples.on)]
            names(samples)[2] <- "subject.id"
            
            samp.pairs <- .samplePairs1(samples)

            ## get original sample filter
            filt.orig <- seqGetFilter(gdsobj)$sample.sel

            var.id <- seqGetData(gdsobj, "variant.id")
            nvar <- length(var.id)
            var.discord <- integer(nvar)
            var.npair <- integer(nvar)
            
            ## for each subject
            nsubj <- nrow(samp.pairs)
            subj.discord <- integer(nsubj)
            subj.nvar <- integer(nsubj)
            for (i in 1:nsubj) {
              if (verbose) message("subject ", i, " of ", nsubj)
              
              ## get genotypes for all samples
              seqSetFilter(gdsobj, sample.id=unlist(samp.pairs[i,], use.names=FALSE),
                           verbose=FALSE)
              geno <- seqGetData(gdsobj, "genotype")
              match <- .genoMatch(geno, check.phase)
              
              ## store concordance by sample in sample list
              subj.discord[i] <- sum(!match, na.rm=TRUE)
              subj.nvar[i] <- sum(!is.na(match))
              
              ## add concordance by variant to variant list
              var.npair <- var.npair + !is.na(match)
              disc <- as.integer(!match)
              disc[is.na(match)] <- 0
              var.discord <- var.discord + disc
            }

            var.df <- data.frame(num.discord=var.discord, num.pair=var.npair,
                                 discord.rate=var.discord/var.npair,
                                 row.names=var.id)

            samp.pairs$num.discord <- subj.discord
            samp.pairs$num.var <- subj.nvar
            samp.pairs$discord.rate <- subj.discord/subj.nvar
            
            ## reset original sample filter
            seqSetFilter(gdsobj, sample.sel=filt.orig, verbose=FALSE)

            list(by.variant=var.df, by.subject=samp.pairs)
            ## })
            }



#####################################################################
# definition for the signature with two SeqVarGDSClass objects
#
#

# returns number of concordant genotypes
.getMatchesConc <- function(geno1, geno2){
  sel <- (geno1 == geno2) & (geno1 != "miss")
  sel
}

# returns index selecting genotype pairs that involve the alt allele
.getAlt <- function(geno1, geno2){
  (geno1 == "het" | geno1 == "alt" | geno2 == "het" | geno2 == "alt") & (geno1 != "miss") & (geno2 != "miss")
}


# returns the number of concordant pairs involving the alt allele
.getMatchesAltConc <- function(geno1, geno2){
  sel <- .getAlt(geno1, geno2) & (geno1 == geno2) & (geno1 != "miss")
  sel
}

.getMatchesHetAlt <- function(geno1, geno2) {
  sel <- (geno1 == "alt" & geno2 == "het") |
    (geno1 == "het" & geno2 == "alt")
  sel
}

.getMatchesHetRef <- function(geno1, geno2) {
  sel <- (geno1 == "ref" & geno2 == "het") |
    (geno1 == "het" & geno2 == "ref")
  sel
}

.getMatchesRefAlt <- function(geno1, geno2){
  sel <- (geno1 == "ref" & geno2 == "alt") | (geno1 == "alt" & geno2 == "ref")
  sel
}


.getMatchesHetHom <- function(geno1, geno2){
  sel <- (geno1 %in% c("ref", "alt") & geno2 %in% c("ref", "alt")) |
    (geno1 == "het" & geno2 == "het")
  sel
}

.getNonMissing <- function(geno1, geno2){
  sel <- (geno1 != "miss") & (geno2 != "miss")
  sel
}


.matchVariants <- function(gr1, gr2, match.on=c("alleles", "position"), allowOverlaps=TRUE){
  
  match.on = match.arg(match.on)
  
  required <- c("variant.id", "ref", "alt", "snv")
  names1 <- names(elementMetadata(gr1))
  if (length(setdiff(required, names1)) > 0) stop(paste("gr1 must have metadata", paste(required, collapse=", ")))
  
  names2 <- names(elementMetadata(gr2))
  if (length(setdiff(required, names2)) > 0) stop(paste("gr2 must have metadata", paste(required, collapse=", ")))
  
  # subset to only biallelic snvs
  #sel.snv <- gr1$snv & gr2$snv
  gr1 <- gr1[gr1$snv]
  gr2 <- gr2[gr2$snv]
  
  # find overlaps
  overlaps <- findOverlaps(gr1, gr2)
  

  
  overlapping.1 <- gr1[queryHits(overlaps)]
  overlapping.2 <- gr2[subjectHits(overlaps)]
  
  # data frame to track overlaps
  overlapping <- data.frame(variant.id.1=overlapping.1$variant.id,
                            variant.id.2=overlapping.2$variant.id,
                            stringsAsFactors=F)

  # check alleles if requested
  sel.same <- (overlapping.1$ref == overlapping.2$ref) & (overlapping.1$alt == overlapping.2$alt)
  sel.flip <- (overlapping.1$ref == overlapping.2$alt) & (overlapping.1$alt == overlapping.2$ref)
  if (match.on == "alleles"){
    overlapping$recode <- sel.flip
    overlapping <- overlapping[sel.same | sel.flip, ]
  } else {
    overlapping$recode <- FALSE
  }
  
  if (!allowOverlaps){
    # choose the first overlapping variant
    overlapping <- overlapping[!duplicated(overlapping$variant.id.1), ]
    overlapping <- overlapping[!duplicated(overlapping$variant.id.2), ]
  }
  
  overlapping
  
}

.matchSamples <- function(samp1, samp2) {
  
  if (!("sample.id" %in% names(samp1)) | !("sample.id" %in% names(samp2))) stop("sample data frames must have sample.id")
  if (!("subject.id" %in% names(samp1)) | !("subject.id" %in% names(samp2))) stop("sample data frames must have subject.id")
  
  # match on subject.ids
  subjects <- intersect(samp1$subject.id, samp2$subject.id)
  
  samp1 <- samp1[samp1$subject.id %in% subjects, ]
  samp2 <- samp2[samp2$subject.id %in% subjects, ]
  
  names(samp1)[names(samp1) == "sample.id"] <- "sample.id.1"
  names(samp2)[names(samp2) == "sample.id"] <- "sample.id.2"
  
  samp <- merge(samp1, samp2)
  samp
  
}


.getGRanges <- function(gds){
  gr <- granges(gds)
  gr$variant.id <- seqGetData(gds, "variant.id")
  gr$ref <- refChar(gds)
  gr$alt <- altChar(gds)
  gr$snv <- isSNV(gds)
  
  gr
}


.getGenotypeClass <- function(x){
  
  # map for genotype classes
  class.map <- c("alt", "het", "ref")
  
  # 0 = alt/alt, 1 = het, 2 = ref/ref, so we can just subset class map by the dosage plus 1
  tmp <- class.map[x + 1]
  tmp[is.na(tmp)] <- "miss"
  
  if (is.matrix(x)) tmp <- matrix(tmp, nrow=nrow(x))
  
  tmp
  
}


# assumes filters are already set, and will reset filters
setMethod("duplicateDiscordance",
          c("SeqVarData", "SeqVarData"),
          function(gdsobj, obj2, match.samples.on=c("subject.id", "subject.id"), match.variants.on=c("alleles", "position"), discordance.type=c("genotype", "hethom"), by.variant=FALSE, verbose=TRUE){
            
            # deal with arguments
            match.variants.on = match.arg(match.variants.on)
            discordance.type = match.arg(discordance.type)
            
            # save original filters
            originalVariants1 <- seqGetData(gdsobj, "variant.id")
            originalSamples1 <- seqGetData(gdsobj, "sample.id")
            
            originalVariants2 <- seqGetData(obj2, "variant.id")
            originalSamples2 <- seqGetData(obj2, "sample.id")
            
            
            # match samples
            # construct sample data frames for matching
            if (verbose) message("matching samples... ", appendLF=FALSE)
            samp1 <- pData(sampleData(gdsobj))[, c("sample.id", match.samples.on[1])]
            names(samp1)[2] <- "subject.id"

            samp2 <- pData(sampleData(obj2))[, c("sample.id", match.samples.on[2])]
            names(samp2)[2] <- "subject.id"

            samples <- .matchSamples(samp1, samp2)
            if (verbose) message(paste(nrow(samples), "pairs identified!"))
            
            # match variants
            if (verbose) message("matching variants... ", appendLF=FALSE)
            gr1 <- .getGRanges(gdsobj)
            gr2 <- .getGRanges(obj2)
            overlappingVariants <- .matchVariants(gr1, gr2, match.on=match.variants.on)
            if (verbose) message(paste(nrow(overlappingVariants), "variant matches identified!"))
            
            # set up results data frame -- can just add columns to the samples data frame
            if (!by.variant){

              nsamp <- nrow(samples)
              n.variants <- integer(nsamp)
              n.concordant <- integer(nsamp)
              
              if (discordance.type == "genotype") {
                n.alt <- integer(nsamp)
                n.alt.conc <- integer(nsamp)
                n.het.ref <- integer(nsamp)
                n.het.alt <- integer(nsamp)
                n.ref.alt <- integer(nsamp)
              }
            } else {

              nvar <- nrow(overlappingVariants)
              n.samples <- integer(nvar)
              n.concordant <- integer(nvar)
              
              if (discordance.type == "genotype"){
                n.alt <- integer(nvar)
                n.alt.conc <- integer(nvar)
                n.het.ref <- integer(nvar)
                n.het.alt <- integer(nvar)
                n.ref.alt <- integer(nvar)
              }
            }
            
            # set filters for the variants -- will still need to order them properly
            seqSetFilter(gdsobj, variant.id=overlappingVariants$variant.id.1, verbose=FALSE)
            seqSetFilter(obj2, variant.id=overlappingVariants$variant.id.2, verbose=FALSE)
            
            for (i in seq_along(samples$subject.id)){
              
              if (verbose & (i %% 10) == 0){
                message(paste("sample pair", i, "of", nsamp))
              }
              
              # set sample filter for this pair
              seqSetFilter(gdsobj, sample.id=samples$sample.id.1[i], verbose=FALSE)
              seqSetFilter(obj2, sample.id=samples$sample.id.2[i], verbose=FALSE)
              
              # prepare genotype data frame
              dos1 <- refDosage(gdsobj)
              dos2 <- refDosage(obj2)
              
              # order genotypes appropriately
              dos1 <- dos1[, as.character(overlappingVariants$variant.id.1)]
              dos2 <- dos2[, as.character(overlappingVariants$variant.id.2)]
              
              # recode if the alt/ref alleles are switched
              dos2[overlappingVariants$recode] <- 2 - dos2[overlappingVariants$recode]
              
              # remove missing genotypes
              sel <- !is.na(dos1) & !is.na(dos2)
              #dos1 <- dos1[sel]
              #dos2 <- dos2[sel]
              
              class1 <- .getGenotypeClass(dos1)
              class2 <- .getGenotypeClass(dos2)
              
              if (!by.variant){
                
                # store information about discordant variants for this pair
                n.variants[i] <- sum(sel)
                
                if (discordance.type == "hethom"){
                  n.concordant[i] <- sum(.getMatchesHetHom(class1, class2)[sel])
                } else {
                  n.concordant[i] <- sum(.getMatchesConc(class1, class2)[sel])
                  n.alt[i] <- sum(.getAlt(class1, class2)[sel])
                  n.alt.conc[i] <- sum(.getMatchesAltConc(class1, class2)[sel])
                  n.het.ref[i] <- sum(.getMatchesHetRef(class1, class2)[sel])
                  n.het.alt[i] <- sum(.getMatchesHetAlt(class1, class2)[sel])
                  n.ref.alt[i] <- sum(.getMatchesRefAlt(class1, class2)[sel])
                }
              } else {
                
                n.samples <- n.samples + sel
                
                if (discordance.type == "hethom"){
                  n.concordant <- n.concordant + .getMatchesHetHom(class1, class2)  
                } else if (discordance.type == "genotype") {
                  n.concordant <- n.concordant + .getMatchesConc(class1, class2)
                  n.alt <- n.alt + .getAlt(class1, class2)
                  n.alt.conc <- n.alt.conc + .getMatchesAltConc(class1, class2)
                  n.het.ref <- n.het.ref + .getMatchesHetRef(class1, class2)
                  n.het.alt <- n.het.alt + .getMatchesHetAlt(class1, class2)
                  n.ref.alt <- n.ref.alt + .getMatchesRefAlt(class1, class2)
                }
              }
              
            }
            
            # reset original filters
            seqSetFilter(gdsobj, sample.id=originalSamples1, variant.id=originalVariants1, verbose=FALSE)
            seqSetFilter(obj2, sample.id=originalSamples2, variant.id=originalVariants2, verbose=FALSE)
            
            # return results data frame
            if (!by.variant){
              samples <- cbind(samples, n.variants, n.concordant)
              if (discordance.type == "genotype") {
                samples <- cbind(samples,
                                 n.alt, n.alt.conc,
                                 n.het.ref, n.het.alt, n.ref.alt)
              }
              return(samples)
            } else {
              overlappingVariants <- cbind(overlappingVariants,
                                           n.samples, n.concordant)
              if (discordance.type == "genotype") {
                overlappingVariants <- cbind(overlappingVariants,
                                             n.alt, n.alt.conc,
                                             n.het.ref, n.het.alt, n.ref.alt)
              }
              return(overlappingVariants)
            }
            
          })


setMethod("duplicateDiscordance",
          c("SeqVarData", "missing"),
          function(gdsobj,  match.samples.on="subject.id", by.variant=FALSE, all.pairs=TRUE, verbose=TRUE){
              
              ## samples should have columns of sample.id, subject.id
              ## find matching sample pairs for subjects (one sample per subject)
              samples <- pData(sampleData(gdsobj))
              if (!(match.samples.on %in% names(samples))) stop(sprintf("%s is not a column in sampleData", match.samples.on))
              samples <- samples[, c("sample.id", match.samples.on)]
              names(samples)[2] <- "subject.id"
              
              samp.pairs <- .samplePairs(samples, all.pairs=all.pairs)

              ## get original sample filter
              filt.orig <- seqGetFilter(gdsobj)$sample.sel

              nsamp <- nrow(samp.pairs)
              
              if (!by.variant) {
                  n.variants <- integer(nsamp)
                  n.concordant <- integer(nsamp)
                  n.alt <- integer(nsamp)
                  n.alt.conc <- integer(nsamp)
                  n.het.ref <- integer(nsamp)
                  n.het.alt <- integer(nsamp)
                  n.ref.alt <- integer(nsamp)
              } else {
                  variant.id <- seqGetData(gdsobj, "variant.id")
                  nvar <- length(variant.id)
                  n.samples <- integer(nvar)
                  n.concordant <- integer(nvar)
                  n.alt <- integer(nvar)
                  n.alt.conc <- integer(nvar)
                  n.het.ref <- integer(nvar)
                  n.het.alt <- integer(nvar)
                  n.ref.alt <- integer(nvar)
              }
              
              for (i in 1:nsamp) {
                  if (verbose) message("sample pair ", i, " of ", nsamp)
                  
                  ## get genotypes for all samples
                  seqSetFilter(gdsobj, sample.id=unlist(samp.pairs[i,-1], use.names=FALSE),
                               verbose=FALSE)
                  geno <- refDosage(gdsobj, use.names=FALSE)
                  dos1 <- geno[1,]
                  dos2 <- geno[2,]
                  rm(geno)
                  
                  sel <- !is.na(dos1) & !is.na(dos2)
                  
                  class1 <- .getGenotypeClass(dos1)
                  class2 <- .getGenotypeClass(dos2)
                  
                  if (!by.variant){
                      # store information about discordant variants for this pair
                      n.variants[i] <- sum(sel)
                      n.concordant[i] <- sum(.getMatchesConc(class1, class2)[sel])
                      n.alt[i] <- sum(.getAlt(class1, class2)[sel])
                      n.alt.conc[i] <- sum(.getMatchesAltConc(class1, class2)[sel])
                      n.het.ref[i] <- sum(.getMatchesHetRef(class1, class2)[sel])
                      n.het.alt[i] <- sum(.getMatchesHetAlt(class1, class2)[sel])
                      n.ref.alt[i] <- sum(.getMatchesRefAlt(class1, class2)[sel])
                  } else {
                      n.samples <- n.samples + sel
                      n.concordant <- n.concordant + .getMatchesConc(class1, class2)
                      n.alt <- n.alt + .getAlt(class1, class2)
                      n.alt.conc <- n.alt.conc + .getMatchesAltConc(class1, class2)
                      n.het.ref <- n.het.ref + .getMatchesHetRef(class1, class2)
                      n.het.alt <- n.het.alt + .getMatchesHetAlt(class1, class2)
                      n.ref.alt <- n.ref.alt + .getMatchesRefAlt(class1, class2)
                  }
              }
              
              ## reset original sample filter
              seqSetFilter(gdsobj, sample.sel=filt.orig, verbose=FALSE)

              # return results data frame
              if (!by.variant){
                  subj.df <- cbind(samp.pairs,
                                   n.variants,
                                   n.concordant,
                                   n.alt,
                                   n.alt.conc,
                                   n.het.ref,
                                   n.het.alt,
                                   n.ref.alt)
                  return(subj.df)
                  
              } else {
                  var.df <- data.frame(variant.id,
                                       n.samples,
                                       n.concordant,
                                       n.alt,
                                       n.alt.conc,
                                       n.het.ref,
                                       n.het.alt,
                                       n.ref.alt)
                  return(var.df)
              }
              
          })



setMethod("duplicateDiscordance",
          c("SeqVarIterator", "missing"),
          function(gdsobj,  match.samples.on="subject.id", by.variant=FALSE, all.pairs=TRUE, verbose=TRUE){
              
              ## samples should have columns of sample.id, subject.id
              ## find matching sample pairs for subjects (one sample per subject)
              samples <- pData(sampleData(gdsobj))
              if (!(match.samples.on %in% names(samples))) stop(sprintf("%s is not a column in sampleData", match.samples.on))
              samples <- samples[, c("sample.id", match.samples.on)]
              names(samples)[2] <- "subject.id"
              
              samp.pairs <- .samplePairs(samples, all.pairs=all.pairs)
              sample.id <- unique(c(samp.pairs$sample.id.1, samp.pairs$sample.id.2))

              ## get original sample filter
              filt.orig <- seqGetFilter(gdsobj)$sample.sel
              seqSetFilter(gdsobj, sample.id=sample.id, verbose=FALSE)
              sample.id <- seqGetData(gdsobj, "sample.id") # in case order is different

              nsamp <- nrow(samp.pairs)
              
              # results
              res <- list()
              n.iter <- length(variantFilter(gdsobj))
              b <- 1
              iterate <- TRUE
              while (iterate) {
                  ## get genotypes for all samples
                  geno <- refDosage(gdsobj, use.names=FALSE)
                  class <- .getGenotypeClass(geno)
                  rownames(class) <- sample.id
                  rm(geno)
                  
                  class1 <- class[samp.pairs$sample.id.1,,drop=FALSE]
                  class2 <- class[samp.pairs$sample.id.2,,drop=FALSE]
                  
                  if (!by.variant){
                      res[[b]] <- cbind(
                          n.variants = rowSums(.getNonMissing(class1, class2)),
                          n.concordant = rowSums(.getMatchesConc(class1, class2)),
                          n.alt = rowSums(.getAlt(class1, class2)),
                          n.alt.conc = rowSums(.getMatchesAltConc(class1, class2)),
                          n.het.ref = rowSums(.getMatchesHetRef(class1, class2)),
                          n.het.alt = rowSums(.getMatchesHetAlt(class1, class2)),
                          n.ref.alt = rowSums(.getMatchesRefAlt(class1, class2))
                      )
                  } else {
                      res[[b]] <- data.frame(
                          variant.id = seqGetData(gdsobj, "variant.id"),
                          n.samples = colSums(.getNonMissing(class1, class2)),
                          n.concordant = colSums(.getMatchesConc(class1, class2)),
                          n.alt = colSums(.getAlt(class1, class2)),
                          n.alt.conc = colSums(.getMatchesAltConc(class1, class2)),
                          n.het.ref = colSums(.getMatchesHetRef(class1, class2)),
                          n.het.alt = colSums(.getMatchesHetAlt(class1, class2)),
                          n.ref.alt = colSums(.getMatchesRefAlt(class1, class2))
                      )
                  }
                  
                  if (verbose & b %% 100 == 0) {
                      message(paste("Iteration", b , "of", n.iter, "completed"))
                  }
                  b <- b + 1
                  iterate <- iterateFilter(gdsobj, verbose=FALSE)
              }
              
              ## reset original sample filter
              seqSetFilter(gdsobj, sample.sel=filt.orig, verbose=FALSE)

              # return results data frame
              if (!by.variant){
                  ## sum over variant blocks
                  res <- Reduce(`+`, res)
                  subj.df <- cbind(samp.pairs, res)
                  return(subj.df)
                  
              } else {
                  var.df <- do.call(rbind, res)
                  return(var.df)
              }
          })
