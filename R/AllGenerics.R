
setGeneric("applyMethod",
           signature=c("gdsobj", "FUN", "variant"),
           function(gdsobj, FUN, variant=NULL, sample=NULL, ...)
             standardGeneric("applyMethod"))


## allele methods
setGeneric("refChar",
           function(gdsobj, ...)
             standardGeneric("refChar"))

setGeneric("altChar",
           function(gdsobj, ...)
             standardGeneric("altChar"))

setGeneric("nAlleles",
           function(gdsobj, ...)
             standardGeneric("nAlleles"))


## descriptive methods
setGeneric("isVariant",
           function(gdsobj, ...)
             standardGeneric("isVariant"))

setGeneric("isSNV",
           function(gdsobj, ...)
             standardGeneric("isSNV"))


## data retrieval and formatting
setGeneric("getGenotype",
           function(gdsobj, ...)
             standardGeneric("getGenotype"))

setGeneric("getGenotypeAlleles",
           function(gdsobj, ...)
             standardGeneric("getGenotypeAlleles"))

setGeneric("refDosage",
           function(gdsobj, ...)
             standardGeneric("refDosage"))

setGeneric("getVariableLengthData",
           function(gdsobj, var.name, ...)
             standardGeneric("getVariableLengthData"))


## metrics
setGeneric("titv",
           function(gdsobj, ...)
             standardGeneric("titv"))

setGeneric("alleleFrequency",
           function(gdsobj, ...)
             standardGeneric("alleleFrequency"))

setGeneric("missingGenotypeRate",
           function(gdsobj, ...)
             standardGeneric("missingGenotypeRate"))

setGeneric("heterozygosity",
           function(gdsobj, ...)
             standardGeneric("heterozygosity"))

setGeneric("homozygosity",
           function(gdsobj, ...)
             standardGeneric("homozygosity"))

setGeneric("meanBySample",
           function(gdsobj, ...)
             standardGeneric("meanBySample"))

## other
setGeneric("pca",
           function(gdsobj, ...)
             standardGeneric("pca"))

setGeneric("hwe",
           function(gdsobj, ...)
             standardGeneric("hwe"))

setGeneric("inbreedCoeff",
           function(gdsobj, ...)
             standardGeneric("inbreedCoeff"))

setGeneric("mendelErr",
           function(gdsobj, ...)
             standardGeneric("mendelErr"))

setGeneric("duplicateDiscordance",
           function(gdsobj, ...)
             standardGeneric("duplicateDiscordance"))
