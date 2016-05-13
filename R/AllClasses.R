setClass("SeqVarData", contains="SeqVarGDSClass",
         slots=c(sampleData = "AnnotatedDataFrame",
                 variantData = "AnnotatedDataFrame"))
