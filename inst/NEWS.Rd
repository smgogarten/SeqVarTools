\name{NEWS}
\title{NEWS for SeqVarTools}
\section{Changes in version 1.5.1}{
  \itemize{
    \item Use existing isSNV generic from VariantAnnotation instead of redefining
    \item Use BiocStyle for vignette
  }
}
\section{Changes in version 1.7.3}{
  \itemize{
    \item Add methods for duplicateDiscordance with two datasets
    \item Add alternateAlleleDetection
  }
}
\section{Changes in version 1.7.4}{
  \itemize{
    \item HWE works on biallelic INDELs as well as SNVs
  }
}
\section{Changes in version 1.7.5}{
  \itemize{
    \item Added SeqVarData class to combine sample annotation with GDS object
  }
}
\section{Changes in version 1.7.6}{
  \itemize{
    \item HWE method returns additional columns and allows permutation of genotypes
  }
}
\section{Changes in version 1.7.7}{
  \itemize{
    \item added methods to calculate and plot reference allele fraction
  }
}
\section{Changes in version 1.7.9}{
  \itemize{
    \item added more options to return dosage of different alleles
  }
}
\section{Changes in version 1.9.2}{
  \itemize{
    \item added by.variant option to duplicateDiscordance for two gds files
  }
}
\section{Changes in version 1.9.4}{
  \itemize{
    \item alleleDosage returns list with dosage of each allele separately
  }
}
\section{Changes in version 1.9.8}{
  \itemize{
    \item duplicateDiscordance and alternateAlleleDetection require SeqVarData objects; both can match on a subject.id instead of sample.id
  }
}
\section{Changes in version 1.9.9}{
  \itemize{
    \item duplicateDiscordance across two GDS files can match on either position or position and alleles, with the ability to recode genotypes if th reference allele in one dataset is the alternate allele in the other dataset
  }
}
\section{Changes in version 1.9.10}{
  \itemize{
    \item duplicateDiscordance across two GDS files can calculate discordance based on heterozygote/homozygote status instead of genotype
  }
}
\section{Changes in version 1.9.11}{
  \itemize{
    \item Add Firth test option to regression
    \item Bug fix for refFracPlot: hets significantly different from 0.5 plotted as triangles, median line shown
  }
}
\section{Changes in version 1.11.1}{
  \itemize{
    \item Add variantData slot to SeqVarData class
  }
}
\section{Changes in version 1.11.3}{
  \itemize{
    \item Add hethom method to calculate heterozygosity / non-reference homozygosity in one step
    \item Add countSingletons method
  }
}
\section{Changes in version 1.15.1}{
  \itemize{
    \item Following SeqArray, remove dependency on VariantAnnotation
    \item Add generic isSNV (replacing previous import of this generic from VariantAnnotation)
  }
}
\section{Changes in version 1.15.2}{
  \itemize{
    \item Added iterator classes: SeqVarBlockIterator, SeqVarRangeIterator, SeqVarWindowIterator, SeqVarListIterator.
    \item Creating a SeqVarData object with missing sample or variant annotation will store 0-column data frames in sampleData or variantData, instead of duplicating sample.id and variant.id.
    \item Added methods to return variant data in expanded form, with one row per alternate allele.
  }
}

\section{Changes in version 1.15.3}{
  \itemize{
    \item alleleFrequency method accounts for sex when computing frequency for X and Y chromosomes.
  }
}

\section{Changes in version 1.17.2}{
  \itemize{
    \item Change implementation of iterator classes to identify indices of selected variants on object creation and store in variantList slot. All iterator classes now extend new class SeqVarIterator.
  }
}
