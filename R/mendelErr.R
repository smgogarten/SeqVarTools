.triosFromPedigree <- function(pedigree) {
  fams <- unique(pedigree$family)
  trios <- list()
  for (f in fams) {
    fam <- pedigree[pedigree$family == f,]
    for (i in 1:nrow(fam)) {
      if (fam$mother[i] %in% fam$individ & fam$father[i] %in% fam$individ) {
        trio <- c(child=fam$sample.id[i],
                  mother=fam$sample.id[fam$individ == fam$mother[i]],
                  father=fam$sample.id[fam$individ == fam$father[i]])
        trios[[length(trios)+1]] <- trio
      } else if (fam$mother[i] %in% fam$individ) {
        trio <- c(child=fam$sample.id[i],
                  mother=fam$sample.id[fam$individ == fam$mother[i]])
        trios[[length(trios)+1]] <- trio
      } else if (fam$father[i] %in% fam$individ) {
        trio <- c(child=fam$sample.id[i],
                  father=fam$sample.id[fam$individ == fam$father[i]])
        trios[[length(trios)+1]] <- trio
      }
    }
  }
  trios
}

.trioType <- function(trio) {
  if ("mother" %in% trio & "father" %in% trio) {
    "both.parents"
  } else if ("mother" %in% trio) {
    "mother.only"
  } else if ("father" %in% trio) {
    "father.only"
  }
}

.alleleMatch <- function(geno, allele, parent) {
  match1 <- geno[allele,"child",] == geno[1,parent,]
  match2 <- geno[allele,"child",] == geno[2,parent,]
  match1 | match2
}

## classes of Mendelian errors
## autosomes, pseudoautosomal region, and X chrom female child
## father  mother  child
## AA      AA      AB
## BB      BB      AB
## BB      **      AA
## **      BB      AA
## AA      **      BB
## **      AA      BB
.autosomeErr <- function(geno) {
  type <- .trioType(dimnames(geno)$sample)
  hom <- geno[1,"child",] == geno[2,"child",]
  if (type == "both.parents") {
    mother1match <- .alleleMatch(geno, 1, "mother")
    mother2match <- .alleleMatch(geno, 2, "mother")
    father1match <- .alleleMatch(geno, 1, "father")
    father2match <- .alleleMatch(geno, 2, "father")
    err.het <- !(mother1match | father1match) | !(mother2match | father2match)
    err.mother <- hom & !(mother1match | mother2match)
    err.father <- hom & !(father1match | father2match)
    err <- err.het | err.mother | err.father
    ## deal with missing values
    if (sum(is.na(err) > 0)) {
      hom <- geno[1,"child",] == geno[2,"child",]
      mother.missing <- is.na(mother1match) | is.na(mother2match)
      err[mother.missing] <- err.father[mother.missing]
      father.missing <- is.na(father1match) | is.na(father2match)
      err[father.missing] <- err.mother[father.missing]
    }
  } else if (type == "mother.only") {
    ## can only test for homozygote errors
    mother1match <- .alleleMatch(geno, 1, "mother")
    mother2match <- .alleleMatch(geno, 2, "mother")
    err <- hom & !(mother1match | mother2match)
  } else if (type == "father.only") {
    ## can only test for homozygote errors
    father1match <- .alleleMatch(geno, 1, "father")
    father2match <- .alleleMatch(geno, 2, "father")
    err <- hom & !(father1match | father2match)
  }
  err[is.na(err)] <- FALSE
  err
}

## X chromosome male child
## father  mother  child
## **      AA      BB
## **      BB      AA
.xMaleErr <- function(geno) {
  type <- .trioType(dimnames(geno)$sample)
  if (type == "both.parents" | type == "mother.only") {
    mother1match <- .alleleMatch(geno, 1, "mother")
    mother2match <- .alleleMatch(geno, 2, "mother")
    err <- !(mother1match | mother2match)
  } else {
    err <- rep(FALSE, dim(geno)[3])
  }
  err[is.na(err)] <- FALSE
  err
}

## Y chromosome male child
## father  mother  child
## AA      **      BB
## BB      **      AA
.yMaleErr <- function(geno) {
  type <- .trioType(dimnames(geno)$sample)
  if (type == "both.parents" | type == "father.only") {
    father1match <- .alleleMatch(geno, 1, "father")
    father2match <- .alleleMatch(geno, 2, "father")
    err <- !(father1match | father2match)
  } else {
    err <- rep(FALSE, dim(geno)[3])
  }
  err[is.na(err)] <- FALSE
  err
}


setMethod("mendelErr",
          "SeqVarGDSClass",
          function(gdsobj, pedigree, use.names=FALSE,
                   autosomes=1:22, xchrom="X", ychrom="Y",
                   verbose=TRUE) {
            ## pedigree shoud have format (family, individ, father, mother, sex, sample.id)
            ## generate trios (and duos) from pedigree
            trios <- .triosFromPedigree(pedigree)
            
            ## deal with duplicate samples (to-do later)
            
            ## chromosome types
            chrom <- seqGetData(gdsobj, "chromosome")
            auto <- chrom %in% autosomes
            x <- chrom %in% xchrom
            y <- chrom %in% ychrom
            ## pseudoautosomal?
            ## mitochondrial?           

            ## get original sample filter
            filt.orig <- seqGetFilter(gdsobj)$sample.sel
            
            ## vector to store errors by variant
            var.err <- integer(.nVar(gdsobj))

            ## for each trio
            trio.err <- integer()
            for (t in trios) {
              ## get child + parent(s) genotype
              seqSetFilter(gdsobj, sample.id=t, verbose=verbose)
              geno <- seqGetData(gdsobj, "genotype")
              sample.id <- seqGetData(gdsobj, "sample.id")
              dimnames(geno)$sample <- names(t)[match(sample.id, t)]
              
              ## add to count of errors 
              child.sex <- pedigree$sex[pedigree$sample.id == t["child"]]
              aut.type <- switch(child.sex, F=(auto | x), M=auto)
              aut.err <- .autosomeErr(geno[,,aut.type])
              var.err[aut.type] <- var.err[aut.type] + aut.err

              if (child.sex == "M") {
                x.err <- .xMaleErr(geno[,,x])
                var.err[x] <- var.err[x] + x.err
                y.err <- .yMaleErr(geno[,,y])
                var.err[y] <- var.err[y] + y.err
              } else {
                x.err <- NULL
                y.err <- NULL
              }

              ## keep track of errors by trio
              trio.err[as.character(t["child"])] <- sum(c(aut.err, x.err, y.err))
            }
              
            ## reset original sample filter
            seqSetFilter(gdsobj, sample.sel=filt.orig, verbose=FALSE)

            if (use.names)
              names(var.err) <- seqGetData(gdsobj, "variant.id")
            list(by.variant=var.err, by.trio=trio.err)
          })

