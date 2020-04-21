## Global variables

#--------- hg19 ref genome ----------#
.onLoad <- function(libname, pkgname){
   DEFAULT_GENOME <- tryCatch({
      BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
   }, error=function(e){ return() })

   if(is.null(DEFAULT_GENOME)){
      warning("
   No reference genome loaded. Please install and load a BSgenome.
   For example:
      install.packages('BiocManager')
      BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
      library('BSgenome.Hsapiens.UCSC.hg19')

   Then specify the BSgenome to the ref.genome arguemnts to the relevant functions.
   For example:
      extractSigsSnv(..., ref.genome=BSgenome.Hsapiens.UCSC.hg19)
")
   }

   assign('DEFAULT_GENOME', DEFAULT_GENOME, envir=parent.env(environment()))
}

#--------- SNV context types ----------#
SUBSTITUTIONS <- c('C>A','C>G','C>T','T>A','T>C','T>G')

C_TRIPLETS <- c(
   "ACA", "ACC", "ACG", "ACT",
   "CCA", "CCC", "CCG", "CCT",
   "GCA", "GCC", "GCG", "GCT",
   "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS <- c(
   "ATA", "ATC", "ATG", "ATT",
   "CTA", "CTC", "CTG", "CTT",
   "GTA", "GTC", "GTG", "GTT",
   "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 <- c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
SUBSTITUTIONS_96 <- rep(SUBSTITUTIONS, each=16)
SUBS_CONTEXTS_96 <- paste0(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3))


