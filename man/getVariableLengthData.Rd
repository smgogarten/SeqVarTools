\name{getVariableLengthData}
\alias{getVariableLengthData}
\alias{getVariableLengthData,SeqVarGDSClass,character-method}

\title{Get variable-length data}

\description{
  Get data with multiple values per sample from a GDS object and return
  as an array
}

\usage{
\S4method{getVariableLengthData}{SeqVarGDSClass,character}(gdsobj, var.name, use.names=TRUE)
}

\arguments{
  \item{gdsobj}{A \code{\link{SeqVarGDSClass}} object with VCF data.}
  \item{var.name}{Character string with name of the variable, most
    likely "annotation/format/VARIABLE_NAME".}
  \item{use.names}{A logical indicating whether to assign sample and
    variant IDs as dimnames of the resulting matrix.}
}

\details{
  Data which are indicated as having variable length (possibly different
  numbers of values for each variant) in the VCF header are stored as
  variable-length data in the GDS file.  Each such data object has two
  components, "length" and "data."  "length" indicates how many values
  there are for each variant, while "data" is a matrix with one row per
  sample and columns defined as all values for variant 1, followed by
  all values for variant 2, etc.

  \code{getVariableLengthData} converts this format to a 3-dimensional
  array, where the length of the first dimension is the maximum number
  of values in "length," and the remaining dimensions are sample and
  variant.  Missing values are given as NA.  If the first dimension of
  this array would have length 1, the result is converted to a matrix.
}

\value{
  An array with dimensions [n, sample, variant] where n is the maximum
  number of values possible for a given sample/variant cell.  If n=1, a
  matrix with dimensions [sample,variant].
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SeqVarGDSClass}},
  \code{\link{applyMethod}},
  \code{\link{seqGetData}}
}

\examples{
file <- system.file("extdata", "gl_chr1.gds", package="SeqVarTools")
gds <- seqOpen(file)
## genotype likelihood 
gl <- seqGetData(gds, "annotation/format/GL")
names(gl)
gl$length
## 3 values per variant - likelihood of RR,RA,AA genotypes
dim(gl$data)
## 85 samples (rows) and 9 variants with 3 values each - 27 columns

gl.array <- getVariableLengthData(gds, "annotation/format/GL")
dim(gl.array)
## 3 genotypes x 85 samples x 9 variants
head(gl.array[1,,])
head(gl.array[2,,])
head(gl.array[3,,])

## genotype dosage
ds <- seqGetData(gds, "annotation/format/DS")
names(ds)
ds$length
## 1 value per variant
dim(ds$data)
## 85 samples (rows) and 9 variants (columns)

ds.array <- getVariableLengthData(gds, "annotation/format/DS")
dim(ds.array)
## 85 samples x 9 variants
head(ds.array)

seqClose(gds)
}

\keyword{manip}
