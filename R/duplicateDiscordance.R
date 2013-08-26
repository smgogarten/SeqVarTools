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
          "SeqVarGDSClass",
          function(gdsobj, samples, check.phase=FALSE, verbose=TRUE) {
            ## samples should have columns of sample.id, subject.id
            ## find matching sample pairs for subjects (one sample per subject)
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
