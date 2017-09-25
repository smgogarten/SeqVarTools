## define pseudoautosomal region

setMethod("chromWithPAR",
          "SeqVarGDSClass",
          function(gdsobj, genome.build=c("hg19", "hg38")) {
              build <- match.arg(genome.build)
              chr <- seqGetData(gdsobj, "chromosome")
              if (build == "hg19") {
                  PAR <- GRanges(seqnames=c("X", "X", "Y", "Y"),
                                 ranges=IRanges(start=c(60001, 154931044, 10001, 59034050),
                                                end=c(2699520, 155260560, 2649520, 59363566)))
              } else if (build == "hg38") {
                  PAR <- GRanges(seqnames=c("X", "X", "Y", "Y"),
                                 ranges=IRanges(start=c(10001, 155701383, 10001, 56887903),
                                                end=c(2781479, 156030895, 2781479, 57217415)))
              }
              chr[.rangesToSel(gdsobj, PAR)] <- "PAR"
              chr
          })
