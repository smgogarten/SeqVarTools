test_hethom_sum <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  nsamp <- SeqVarTools:::.nSamp(gds)
  nvar <- SeqVarTools:::.nVar(gds)
  
  miss <- missingGenotypeRate(gds)
  n.miss <- miss*nsamp
  n.nonmiss <- (1-miss)*nsamp
  checkEquals(rep(nsamp, nvar),
              heterozygosity(gds)*n.nonmiss + homozygosity(gds)*n.nonmiss + n.miss)
  
  miss <- missingGenotypeRate(gds, margin="by.sample")
  n.miss <- miss*nvar
  n.nonmiss <- (1-miss)*nvar
  checkEquals(rep(nvar, nsamp),
              heterozygosity(gds, margin="by.sample")*n.nonmiss +
              homozygosity(gds, margin="by.sample")*n.nonmiss +
              n.miss)
  seqClose(gds)
}

test_hom_sum <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  checkEquals(homozygosity(gds, allele="any"),
              homozygosity(gds, allele="ref") + homozygosity(gds, allele="alt"))
  checkEquals(homozygosity(gds, allele="any", margin="by.sample"),
              homozygosity(gds, allele="ref", margin="by.sample") +
              homozygosity(gds, allele="alt", margin="by.sample"))
  seqClose(gds)
}

test_heterozygosity <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  checkIdentical(colSums(geno[1,,] != geno[2,,], na.rm=TRUE) /
                 colSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 heterozygosity(gds, margin="by.variant"))
  checkIdentical(rowSums(geno[1,,] != geno[2,,], na.rm=TRUE) /
                 rowSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 heterozygosity(gds, margin="by.sample"))
  seqClose(gds)
}

test_homozygosity_any <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  hom <- geno[1,,] == geno[2,,]
  checkIdentical(colSums(hom, na.rm=TRUE) /
                 colSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="any", margin="by.variant"))
  checkIdentical(rowSums(hom, na.rm=TRUE) /
                 rowSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="any", margin="by.sample"))
  seqClose(gds)
}

test_homozygosity_ref <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  hom <- geno[1,,] == geno[2,,] & geno[1,,] == 0
  checkIdentical(colSums(hom, na.rm=TRUE) /
                 colSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="ref", margin="by.variant"))
  checkIdentical(rowSums(hom, na.rm=TRUE) /
                 rowSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="ref", margin="by.sample"))
  seqClose(gds)
}

test_homozygosity_alt <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  geno <- seqGetData(gds, "genotype")
  hom <- geno[1,,] == geno[2,,] & geno[1,,] > 0
  checkIdentical(colSums(hom, na.rm=TRUE) /
                 colSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="alt", margin="by.variant"))
  checkIdentical(rowSums(hom, na.rm=TRUE) /
                 rowSums(!is.na(geno[1,,]) & !is.na(geno[2,,])),
                 homozygosity(gds, allele="alt", margin="by.sample"))
  seqClose(gds)
}

test_hethom <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  hh <- hethom(gds)
  geno <- seqGetData(gds, "genotype")
  hom <- geno[1,,] == geno[2,,] & geno[1,,] > 0
  checkIdentical(rowSums(geno[1,,] != geno[2,,], na.rm=TRUE) /
                 rowSums(hom, na.rm=TRUE),
                 hh)
  checkEquals(heterozygosity(gds, margin="by.sample") /
              homozygosity(gds, allele="alt", margin="by.sample"),
              hh)
  seqClose(gds)
}

test_heterozygosity_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  het.var <- heterozygosity(gds, margin="by.variant")
  het.samp <- heterozygosity(gds, margin="by.sample")
  seqSetFilter(gds)
  checkIdentical(het.var,
                 applyMethod(gds, heterozygosity, variant=var.id, sample=samp.id,
                             margin="by.variant"))
  checkIdentical(het.samp,
                 applyMethod(gds, heterozygosity, variant=var.id, sample=samp.id,
                             margin="by.sample"))
  seqClose(gds)
}

test_homozygosity_apply <- function() {
  gds <- seqOpen(seqExampleFileName("gds"))
  var.id <- 101:110
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  hom.var <- homozygosity(gds, margin="by.variant")
  hom.samp <- homozygosity(gds, margin="by.sample")
  seqSetFilter(gds)
  checkIdentical(hom.var,
                 applyMethod(gds, homozygosity, variant=var.id, sample=samp.id,
                             margin="by.variant"))
  checkIdentical(hom.samp,
                 applyMethod(gds, homozygosity, variant=var.id, sample=samp.id,
                             margin="by.sample"))
  seqClose(gds)
}
