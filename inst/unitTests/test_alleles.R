test_parseRefAllele <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T")
  checkIdentical(c("A","AA","A","A"), SeqVarTools:::.parseRefAllele(x))
}

test_parseAltAllele_n0 <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T")
  checkIdentical(c("G","G","GG","G,T"), SeqVarTools:::.parseAltAllele(x, n=0))
}

test_parseAltAllele_n1 <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T", "A,G,T,C", "A,G,TT,C")
  checkIdentical(c("G","G","GG","G","G","G"), SeqVarTools:::.parseAltAllele(x, n=1))
}

test_parseAltAllele_n2 <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T", "A,G,T,C", "A,G,TT,C")
  checkIdentical(c(NA,NA,NA,"T","T","TT"), SeqVarTools:::.parseAltAllele(x, n=2))
}

test_parseAltAllele_n3 <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T", "A,G,T,C", "A,G,TT,C")
  checkIdentical(c(NA,NA,NA,NA,"C","C"), SeqVarTools:::.parseAltAllele(x, n=3))
}

## test_parseNumAlleles <- function() {
##   x <- c("A", "A,A", "A,A,A", "A,A,A,A")
##   checkIdentical(1:4, SeqVarTools:::.parseNumAlleles(x))
## }

test_maxAlleleLength <- function() {
  x <- c("A,G", "AA,G", "A,GG", "A,G,T")
  checkIdentical(as.integer(c(1,2,2,1)), SeqVarTools:::.maxAlleleLength(x))
}
