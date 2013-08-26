test_parseVariableLength <- function() {
  x <- list(length=c(3,3,2,1,3),
            data=matrix(c(1,2,3,1,2,3,1,2,1,1,2,3), nrow=4, ncol=12, byrow=TRUE))
  y <- SeqVarTools:::.parseVariableLength(x)
  for (n in 1:3) checkTrue(all(y[n,,] %in% c(n,NA)))
  checkTrue(all(is.na(y[2,,4])))
  checkTrue(all(is.na(y[3,,3:4])))
}

test_parseVariableLength_length1 <- function() {
  x <- list(length=c(1,1,1,1,1),
            data=matrix(1:5, nrow=4, ncol=5, byrow=TRUE))
  checkIdentical(x$data, SeqVarTools:::.parseVariableLength(x))
}

test_parseVariableLength_length0 <- function() {
  x <- list(length=c(2,2,2,0,2),
            data=matrix(rep(1:2,4), nrow=3, ncol=8, byrow=TRUE))
  y <- array(rep(1:2,15), dim=c(2,3,5))
  y[,,4] <- NA
  checkIdentical(y, SeqVarTools:::.parseVariableLength(x))
}

test_parseVariableLength_length01 <- function() {
  x <- list(length=c(1,1,1,0,1),
            data=matrix(1:4, nrow=3, ncol=4, byrow=TRUE))
  y <- matrix(c(1,2,3,NA,4), nrow=3, ncol=5, byrow=TRUE)
  checkEquals(y, SeqVarTools:::.parseVariableLength(x))
}

test_getVariableLength_apply <- function() {
  gds <- seqOpen(system.file("extdata", "gl_chr1.gds", package="SeqVarTools"))
  var.id <- 2:6
  samp.id <- seqGetData(gds, "sample.id")[6:10]
  seqSetFilter(gds, variant.id=var.id, sample.id=samp.id)
  vl <- getVariableLengthData(gds, "annotation/format/GL")
  seqSetFilter(gds)
  checkIdentical(vl,
                 applyMethod(gds, getVariableLengthData, variant=var.id,
                             sample=samp.id, var.name="annotation/format/GL"))
  seqClose(gds)
}
