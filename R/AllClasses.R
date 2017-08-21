setClass("SeqVarData", contains="SeqVarGDSClass",
         slots=c(sampleData = "AnnotatedDataFrame",
                 variantData = "AnnotatedDataFrame"))

setClass("SeqVarBlockIterator", contains="SeqVarData",
         slots=c(variantBlock="integer", # number of variants
                 lastVariant="environment")) # allow pass-by-reference for this slot

setClass("SeqVarRangeIterator", contains="SeqVarData",
         slots=c(variantRanges="GRanges",
                 lastRange="environment"))

setClass("SeqVarWindowIterator", contains="SeqVarRangeIterator",
         slots=c(windowSize="integer", # base pairs
                 windowShift="integer"))
