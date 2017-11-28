setClass("SeqVarData", contains="SeqVarGDSClass",
         slots=c(sampleData = "AnnotatedDataFrame",
                 variantData = "AnnotatedDataFrame"))

setClass("SeqVarIterator", contains="SeqVarData",
         slots=c(variantFilter="list",
                 lastFilter="environment")) # allow pass-by-reference for this slot

setClass("SeqVarBlockIterator", contains="SeqVarIterator",
         slots=c(variantBlock="integer")) # number of variants

setClass("SeqVarRangeIterator", contains="SeqVarIterator",
         slots=c(variantRanges="GRanges"))

setClass("SeqVarWindowIterator", contains="SeqVarRangeIterator",
         slots=c(windowSize="integer", # base pairs
                 windowShift="integer"))

setClass("SeqVarListIterator", contains="SeqVarIterator",
         slots=c(variantRanges="GRangesList"))
