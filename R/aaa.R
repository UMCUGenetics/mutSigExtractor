## Global variables


.onLoad <- function(libname, pkgname){

   ## hg19 ref genome --------------------------------
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


## SNV context types --------------------------------
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


## Indel context types --------------------------------
INDEL_CONTEXTS <- (function(){
   rep_len_del <- c(1:5,'6+')
   rep_len_ins <- c(0:4,'5+')
   big_indel_len <- c(2:4,'5+')

   homopolymer_contexts <- c(
      paste0('del.1.C.',rep_len_del),
      paste0('del.1.T.',rep_len_del),
      paste0('ins.1.C.',rep_len_ins),
      paste0('ins.1.T.',rep_len_ins)
   )

   del_rep_contexts <- paste0(
      'del.',rep(big_indel_len,each=length(rep_len_del)),'.rep.',rep(rep_len_del,length(big_indel_len))
   )
   ins_rep_contexts <- paste0(
      'ins.',rep(big_indel_len,each=length(rep_len_del)),'.rep.',rep(rep_len_ins,length(big_indel_len))
   )

   mh_len <- c(1,1:2,1:3, 1:4,'5+')
   mh_del_len <- c(2,3,3,4,4,4,'5+','5+','5+','5+','5+')
   del_mh_contexts <- paste0('del.',mh_del_len,'.mh.',mh_len)

   c(homopolymer_contexts,del_rep_contexts,ins_rep_contexts,del_mh_contexts)
})()
