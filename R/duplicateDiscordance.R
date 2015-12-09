.samplePairs <- function(samples) {
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

setMethod("duplicateDiscordance",
          c("SeqVarData", "missing"),
        function(gdsobj, match.samples.on="subject.id", check.phase=FALSE, verbose=TRUE) {
            ## samples should have columns of sample.id, subject.id
            ## find matching sample pairs for subjects (one sample per subject)
            samples <- pData(sampleData(gdsobj))
            if (!(match.samples.on %in% names(samples))) stop(sprintf("%s is not a column in sampleData", match.samples.on))
            samples <- samples[, c("sample.id", match.samples.on)]
            names(samples)[2] <- "subject.id"
            
            samp.pairs <- .samplePairs(samples)

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
            seqSetFilter(gdsobj, samp.sel=filt.orig, verbose=FALSE)

            list(by.variant=var.df, by.subject=samp.pairs)
          })



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

              samples$n.variants <- NA
              samples$n.concordant <- NA
              
              if (discordance.type == "genotype") {
                samples$n.alt <- NA
                samples$n.alt.conc <- NA
                samples$n.het.ref <- NA
                samples$n.het.alt <- NA
                samples$n.ref.alt <- NA
              }
            } else {
              
              overlappingVariants$n.samples <- 0
              overlappingVariants$n.concordant <- 0
              
              if (discordance.type == "genotype"){
                overlappingVariants$n.alt <- 0
                overlappingVariants$n.alt.conc <- 0
                overlappingVariants$n.het.ref <- 0
                overlappingVariants$n.het.alt <- 0
                overlappingVariants$n.ref.alt <- 0
              }
            }
            
            # set filters for the variants -- will still need to order them properly
            seqSetFilter(gdsobj, variant.id=overlappingVariants$variant.id.1, verbose=FALSE)
            seqSetFilter(obj2, variant.id=overlappingVariants$variant.id.2, verbose=FALSE)
            
            for (i in seq_along(samples$subject.id)){
              
              if (verbose & (i %% 10) == 0){
                message(paste("sample pair", i, "out of", nrow(samples)))
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
                samples$n.variants[i] <- sum(sel)
                
                if (discordance.type == "hethom"){
                  samples$n.concordant[i] <- sum(.getMatchesHetHom(class1, class2)[sel])
                } else {
                  samples$n.concordant[i] <- sum(.getMatchesConc(class1, class2)[sel])
                  samples$n.alt[i] <- sum(.getAlt(class1, class2)[sel])
                  samples$n.alt.conc[i] <- sum(.getMatchesAltConc(class1, class2)[sel])
                  samples$n.het.ref[i] <- sum(.getMatchesHetRef(class1, class2)[sel])
                  samples$n.het.alt[i] <- sum(.getMatchesHetAlt(class1, class2)[sel])
                  samples$n.ref.alt[i] <- sum(.getMatchesRefAlt(class1, class2)[sel])
                }
              } else {
                
                overlappingVariants$n.samples <- overlappingVariants$n.samples + sel
                
                if (discordance.type == "hethom"){
                  overlappingVariants$n.concordant <- overlappingVariants$n.concordant + .getMatchesHetHom(class1, class2)  
                } else if (discordance.type == "genotype") {
                  overlappingVariants$n.concordant <- overlappingVariants$n.concordant + .getMatchesConc(class1, class2)
                  overlappingVariants$n.alt <- overlappingVariants$n.alt + .getAlt(class1, class2)
                  overlappingVariants$n.alt.conc <- overlappingVariants$n.alt.conc + .getMatchesAltConc(class1, class2)
                  overlappingVariants$n.het.ref <- overlappingVariants$n.het.ref + .getMatchesHetRef(class1, class2)
                  overlappingVariants$n.het.alt <- overlappingVariants$n.het.alt + .getMatchesHetAlt(class1, class2)
                  overlappingVariants$n.ref.alt <- overlappingVariants$n.ref.alt + .getMatchesRefAlt(class1, class2)
                }
              }
              
            }
            
            # reset original filters
            seqSetFilter(gdsobj, sample.id=originalSamples1, variant.id=originalVariants1, verbose=FALSE)
            seqSetFilter(obj2, sample.id=originalSamples2, variant.id=originalVariants2, verbose=FALSE)
            
            # return results data frame
            if (by.variant){
              return(overlappingVariants)
            } else {
              return(samples)
            }
            
          })