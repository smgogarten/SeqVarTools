.parseVariableLength <- function(x) {
  if (all(x$length == 1)) {
    x$data
  } else {
    var <- array(dim=c(max(x$length), nrow(x$data), length(x$length)))

    ## assign each element of length to an index of first array dimension
    n.ind <- rep(NA, ncol(x$data))
    j <- 1
    for (i in 1:length(x$length)) {
      len <- x$length[i]
      if (len > 0) {
        n.ind[j:(j + len - 1)] <- 1:len
        j <- j + len
      }
    }

    ## for each index of first array dimension, get values
    for (n in 1:dim(var)[1]) {
      var.ind <- which(x$length >= n)
      var[n,,var.ind] <- x$data[,n.ind == n]
    }

    ## if first array dimension is 1, simplify to a matrix
    if (dim(var)[1] == 1) {
      var <- var[1,,]
    }
    var
  }
}

test_parseVariableLength <- function() {
  x <- list(length=c(3,3,2,1,3),
            data=matrix(c(1,2,3,1,2,3,1,2,1,1,2,3), nrow=4, ncol=12, byrow=TRUE))
  y <- .parseVariableLength(x)
  for (n in 1:3) checkTrue(all(y[n,,] %in% c(n,NA)))
  checkTrue(all(is.na(y[2,,4])))
  checkTrue(all(is.na(y[3,,3:4])))
}

test_parseVariableLength_length1 <- function() {
  x <- list(length=c(1,1,1,1,1),
            data=matrix(1:5, nrow=4, ncol=5, byrow=TRUE))
  checkIdentical(x$data, .parseVariableLength(x))
}

test_parseVariableLength_length0 <- function() {
  x <- list(length=c(2,2,2,0,2),
            data=matrix(rep(1:2,4), nrow=3, ncol=8, byrow=TRUE))
  y <- array(rep(1:2,15), dim=c(2,3,5))
  y[,,4] <- NA
  checkIdentical(y, .parseVariableLength(x))
}

test_parseVariableLength_length01 <- function() {
  x <- list(length=c(1,1,1,0,1),
            data=matrix(1:4, nrow=3, ncol=4, byrow=TRUE))
  y <- matrix(c(1,2,3,NA,4), nrow=3, ncol=5, byrow=TRUE)
  checkEquals(y, .parseVariableLength(x))
}

test_getVariableLength_DP <- function() {
  gds <- SeqVarTools:::.testData()
  checkIdentical(.parseVariableLength(seqGetData(gds, "annotation/format/DP")),
                 getVariableLengthData(gds, "annotation/format/DP", use.names=FALSE))
  seqClose(gds)
}

test_getVariableLength_GL <- function() {
  gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
  gds <- seqOpen(system.file("extdata", "gl_chr1.gds", package="SeqVarTools"))
  gl <- getVariableLengthData(gds, "annotation/format/GL", use.names=FALSE)
  dimnames(gl) <- NULL
  checkIdentical(.parseVariableLength(seqGetData(gds, "annotation/format/GL")), gl)
  seqClose(gds)
}

test_getVariableLength_apply <- function() {
  gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
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

test_getVariableLength_AD <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    ad <- getVariableLengthData(gds, "annotation/format/AD", use.names=FALSE)
    dimnames(ad) <- NULL
    checkIdentical(.parseVariableLength(seqGetData(gds, "annotation/format/AD")), ad)
    seqClose(gds)
}

test_getVariableLength_AB <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "hapmap_exome_chr22.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    ab <- getVariableLengthData(gds, "annotation/format/AB", use.names=FALSE)
    dimnames(ab) <- NULL
    checkIdentical(.parseVariableLength(seqGetData(gds, "annotation/format/AB")), ab)
    seqClose(gds)
}
