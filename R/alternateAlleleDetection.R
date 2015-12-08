# this function uses a lot of subfunctions defined in duplicateDiscordance.R

# some additional subfunctions need to be defined here.
# geno1 = sequence
# geno2 = array
.truePos <- function(geno1, geno2) {
  2*(geno1 == "alt" & geno2 == "alt") +
    (geno1 == "alt" & geno2 == "het") +
    (geno1 == "het" & geno2 == "alt") +
    (geno1 == "het" & geno2 == "het")
}

.trueNeg <- function(geno1, geno2) {
  (geno1 == "het" & geno2 == "het") +
    (geno1 == "het" & geno2 == "ref") +
    (geno1 == "ref" & geno2 == "het") +
    2*(geno1 == "ref" & geno2 == "ref")
}

.falsePos <- function(geno1, geno2) {
  (geno1 == "alt" & geno2 == "het") +
    2*(geno1 == "alt" & geno2 == "ref") +
    (geno1 == "het" & geno2 == "ref")
}

.falseNeg <- function(geno1, geno2) {
  
  (geno1 == "het" & geno2 == "alt") +
    2*(geno1 == "ref" & geno2 == "alt") +
    (geno1 == "ref" & geno2 == "het")
}



# assumes filters are already set, and will reset filters
setMethod("alternateAlleleDetection",
          c("SeqVarData", "SeqVarData"),
          function(gdsobj, gdsobj2, match.samples.on=c("subject.id", "subject.id"), verbose=TRUE){
            
            # save original filters
            originalVariants1 <- seqGetData(gdsobj, "variant.id")
            originalSamples1 <- seqGetData(gdsobj, "sample.id")
            
            originalVariants2 <- seqGetData(gdsobj2, "variant.id")
            originalSamples2 <- seqGetData(gdsobj2, "sample.id")
            
            
            # match samples -- no subject matching for now
            if (verbose) message("matching samples... ", appendLF=FALSE)
            samp1 <- pData(sampleData(gdsobj))[, c("sample.id", match.samples.on[1])]
            names(samp1)[2] <- "subject.id"
            
            samp2 <- pData(sampleData(gdsobj2))[, c("sample.id", match.samples.on[2])]
            names(samp2)[2] <- "subject.id"
            
            samples <- .matchSamples(samp1, samp2)
            if (verbose) message(paste(nrow(samples), "pairs identified!"))
            
            # match variants
            if (verbose) message("matching variants... ", appendLF=FALSE)
            gr1 <- .getGRanges(gdsobj)
            gr2 <- .getGRanges(gdsobj2)
            
            overlappingVariants <- .matchVariants(gr1, gr2, allowOverlaps=FALSE)
            
            if (verbose) message(paste(nrow(overlappingVariants), "non-overlapping variant matches identified!"))
            
            # set filters for the variants -- will still need to order them properly
            seqSetFilter(gdsobj, variant.id=overlappingVariants$variant.id.1, verbose=FALSE)
            seqSetFilter(gdsobj2, variant.id=overlappingVariants$variant.id.2, verbose=FALSE)
            
            overlappingVariants$n.samples <- 0
            overlappingVariants$true.pos <- 0
            overlappingVariants$true.neg <- 0
            overlappingVariants$false.pos <- 0
            overlappingVariants$false.neg <- 0
            
            for (i in seq_along(samples$subject.id)){
              
              if (verbose & (i %% 10) == 0){
                message(paste("sample pair", i, "out of", nrow(samples)))
              }
              
              # set sample filter for this pair
              seqSetFilter(gdsobj, sample.id=samples$sample.id.1[i], verbose=FALSE)
              seqSetFilter(gdsobj2, sample.id=samples$sample.id.2[i], verbose=FALSE)
              
              # prepare genotype data frame
              dos1 <- refDosage(gdsobj)
              dos2 <- refDosage(gdsobj2)
              
              # order genotypes appropriately
              dos1 <- dos1[, as.character(overlappingVariants$variant.id.1)]
              dos2 <- dos2[, as.character(overlappingVariants$variant.id.2)]
              
              # recoded variants where ref in one = alt in other
              dos2[overlappingVariants$recode] <- 2 - dos2[overlappingVariants$recode]
              
              ## remove missing genotypes
              #sel <- !is.na(dos1) & !is.na(dos2)
              #dos1 <- dos1[sel]
              #dos2 <- dos2[sel]
              
              class1 <- .getGenotypeClass(dos1)
              class2 <- .getGenotypeClass(dos2)
              
              nonmiss <- (!is.na(dos1) & !is.na(dos2))
              overlappingVariants$n.samples <- overlappingVariants$n.samples + nonmiss
              overlappingVariants$true.pos <- overlappingVariants$true.pos + .truePos(class1, class2)
              overlappingVariants$true.neg <- overlappingVariants$true.neg + .trueNeg(class1, class2)
              overlappingVariants$false.pos <- overlappingVariants$false.pos + .falsePos(class1, class2)
              overlappingVariants$false.neg <- overlappingVariants$false.neg + .falseNeg(class1, class2)
              
            }
            
            # reset original filters
            seqSetFilter(gdsobj, sample.id=originalSamples1, variant.id=originalVariants1, verbose=FALSE)
            seqSetFilter(gdsobj2, sample.id=originalSamples2, variant.id=originalVariants2, verbose=FALSE)
            
            # return results data frame
            overlappingVariants
            
          })