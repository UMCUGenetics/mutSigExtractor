#' Extract signatures in the format compatible with the random forest HRD classifier
#'
#' @param vcf.snv Path to the vcf file containing SNVs
#' @param vcf.indel Path to the vcf file containing indels
#' @param vcf.sv Path to the vcf file containing SVs (Manta vcf)
#' @param sample.name The name of the sample as a character. Defaults to 'sample' if none is provided.
#'
#' @return A 1-row data frame containing the mutational signature contributions
#' @export
#'
#' @examples
#' library(randomForest)
#' sigs <- extractSigsForHrdClassifier(
#'    vcf.snv = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010/analysis/Arne/DR-10_update/scripts/SNV2//XXXXXXXX.vcf.gz',
#'    vcf.indel = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010/analysis/Arne/DR-10_update/scripts/INDEL2//XXXXXXXX.vcf.gz',
#'    vcf.sv = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data//170217_HMFregXXXXXXXX/XXXXXXXX.vcf.gz',
#'    sample.name = 'XXXXXXXX'
#' )
#' rf_model <- readRDS('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mltoolkit/example/data/rf_hrd_predict.rds')
#' predict(rf_model, sigs, type = 'prob')

extractSigsForHrdClassifier <- function(vcf.snv, vcf.indel, vcf.sv, sample.name = 'sample', sv.caller = 'manta'){

   sigs_snv <- extractSigsSnv(vcf.snv, vcf.filter = 'PASS')
   sigs_indel <- extractSigsIndel(vcf.indel, vcf.filter = 'PASS')
   sigs_sv <- extractSigsSv(vcf.sv, vcf.filter = 'PASS', sv.caller = sv.caller)

   #--------- SNV / SV signatures post-processing ---------#
   ## Rename snv and sv signatures
   rownames(sigs_snv) <- paste0('e.',1:30)
   rownames(sigs_sv) <- paste0('SV',1:6)

   #--------- Indel signatures post-processing ---------#
   ## Flatten length (len) / bases in microhomology (bimh) bins
   indel_sig_contexts <- c('del.rep', 'ins.rep','del.mh','ins.mh','del.none','ins.none')
   sigs_indel_merged <- lapply(indel_sig_contexts,function(i){
      sum(
         sigs_indel[grep(i, rownames(sigs_indel)),]
      )
   })
   sigs_indel_merged <- unlist(sigs_indel_merged)
   names(sigs_indel_merged) <- indel_sig_contexts

   ## Merge del.none and ins.none
   sigs_indel_merged <- c(
      sigs_indel_merged[c('del.rep', 'ins.rep','del.mh','ins.mh')],
      indel.none = sum(sigs_indel_merged[c('del.none','ins.none')])
   )

   ## Convert absolute indel sig contribution to relative indel sig contribution
   sigs_indel_merged <- sigs_indel_merged/sum(sigs_indel_merged)
   sigs_indel_merged <- as.data.frame(sigs_indel_merged)

   #--------- Make final matrix ---------#
   ## Need to make colnames the same, otherwise rbind() will give an error
   colnames(sigs_snv) <- sample.name
   colnames(sigs_indel_merged) <- sample.name
   colnames(sigs_sv) <- sample.name

   return(
      t(rbind(sigs_snv, sigs_indel_merged, sigs_sv))
   )
}
