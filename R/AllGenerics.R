
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

setGeneric("chromWithPAR",
           function(gdsobj, ...)
             standardGeneric("chromWithPAR"))


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

setGeneric("altDosage",
           function(gdsobj, ...)
             standardGeneric("altDosage"))

setGeneric("alleleDosage",
           function(gdsobj, n, ...)
             standardGeneric("alleleDosage"))

setGeneric("expandedAltDosage",
           function(gdsobj, ...)
             standardGeneric("expandedAltDosage"))

setGeneric("expandedVariantIndex",
           function(gdsobj, ...)
             standardGeneric("expandedVariantIndex"))

setGeneric("variantInfo",
           function(gdsobj, ...)
             standardGeneric("variantInfo"))

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

setGeneric("hethom",
           function(gdsobj, ...)
             standardGeneric("hethom"))

setGeneric("meanBySample",
           function(gdsobj, ...)
             standardGeneric("meanBySample"))

setGeneric("countSingletons",
           function(gdsobj, ...)
             standardGeneric("countSingletons"))

setGeneric("refFrac",
           function(gdsobj, ...)
             standardGeneric("refFrac"))

setGeneric("refFracOverHets",
           function(gdsobj, ...)
             standardGeneric("refFracOverHets"))

setGeneric("refFracPlot",
           function(gdsobj, ...)
             standardGeneric("refFracPlot"))


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
           function(gdsobj, obj2, ...)
             standardGeneric("duplicateDiscordance"))

setGeneric("alternateAlleleDetection",
           function(gdsobj, gdsobj2, ...)
             standardGeneric("alternateAlleleDetection"))

setGeneric("regression",
           function(gdsobj, ...)
             standardGeneric("regression"))


## SeqVarData
setGeneric("sampleData",
           function(x)
             standardGeneric("sampleData"))

setGeneric("sampleData<-",
           function(x, value)
             standardGeneric("sampleData<-"))

setGeneric("variantData",
           function(x)
             standardGeneric("variantData"))

setGeneric("variantData<-",
           function(x, value)
             standardGeneric("variantData<-"))

setGeneric("validateSex",
           function(x)
             standardGeneric("validateSex"))


## Iterators
setGeneric("variantBlock",
           function(x)
               standardGeneric("variantBlock"))

setGeneric("lastVariant",
           function(x)
               standardGeneric("lastVariant"))

setGeneric("lastVariant<-",
           function(x, value)
               standardGeneric("lastVariant<-"))

setGeneric("variantRanges",
           function(x)
               standardGeneric("variantRanges"))

setGeneric("lastRange",
           function(x)
               standardGeneric("lastRange"))

setGeneric("lastRange<-",
           function(x, value)
               standardGeneric("lastRange<-"))

setGeneric("iterateFilter",
           function(x, ...)
               standardGeneric("iterateFilter"))

setGeneric("restoreFilter",
           function(x, ...)
               standardGeneric("restoreFilter"))
