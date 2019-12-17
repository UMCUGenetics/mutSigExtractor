## Global variables

#--------- Auto-load hg19 ref genome if it exists ----------#
DEFAULT_GENOME <- 'BSgenome.Hsapiens.UCSC.hg19'

.onLoad <- function(libname, pkgname){
   if(DEFAULT_GENOME %in% rownames(installed.packages())){
      message('Loading the default reference genome: ', DEFAULT_GENOME)
      do.call('library', list(DEFAULT_GENOME))
   } else {
      warning("
         No reference genome loaded. Please install and/or load a BSgenome. For example:

         install.packages('BiocManager')
         BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
         library('BSgenome.Hsapiens.UCSC.hg19')
      ")
   }
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


