#' Simplify and normalize contexts
#'
#' @description Exprimental; use at your own risk!
#'
#' @param contexts A list containing data frames of snv, indel, and/or sv contexts. Alternatively,
#' These can be specified with the arguments: snv, indel, or sv
#' @param snv See argument 'contexts'.
#' @param indel See argument 'contexts'.
#' @param sv See argument 'contexts'.
#' @param simplify.types Which types to flatten. Accepts: 'snv','indel','sv'
#' @param rel.types Which types to convert to relative contribution. Accepts: 'snv','indel','sv'
#'
#' @return A matrix or data frame
#' @export
#'
#' @examples
#' base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Breast_Organoids/matrices/'
#' contexts <- list(
#'    snv = readSigsAsDf(paste0(base_dir,'/snv_contexts')),
#'    indel = readSigsAsDf(paste0(base_dir,'/indel')),
#'    sv = readSigsAsDf(paste0(base_dir,'/sv_contexts'))
#' )
#' transformContexts(contexts, simplify.types = c('snv','indel','sv'), rel.types = c('snv','indel','sv'))

transformContexts <- function(contexts=NULL, snv=NULL, indel=NULL, sv=NULL,
                              simplify.types=NULL, lsqnonneg.types=NULL, rel.types=NULL){

   if(!is.null(contexts)){
      snv <- if(!is.null(contexts$snv)){ contexts$snv }
      indel <- if(!is.null(contexts$indel)){ contexts$indel }
      sv <- if(!is.null(contexts$sv)){ contexts$sv }
   }

   out <- list()

   if(!is.null(snv)){
      if('snv' %in% simplify.types){
         snv_split <- splitDfRegex(SUBSTITUTIONS, df = snv)
         snv <- do.call(cbind, lapply(snv_split, rowSums))
         colnames(snv) <- str_replace_all(SUBSTITUTIONS,'>','.')
      }

      if('snv' %in% lsqnonneg.types){
         snv <- t(apply(contexts$snv, 1, function(i){
            fitToSignatures(SNV_SIGNATURE_PROFILES, i)$x
         }))
         colnames(snv) <- paste0('e.',1:30)
      }

      if('snv' %in% rel.types){
         snv <- snv/rowSums(snv)
      }

      out[['snv']] <- snv
   }

   if(!is.null(indel)){
      if('indel' %in% simplify.types){
         indel_types <- c('del.rep', 'ins.rep', 'del.mh', 'ins.mh', 'del.none', 'ins.none')
         indel_split <- splitDfRegex(indel_types, df = indel)
         indel <- do.call(cbind, lapply(indel_split, rowSums))
         colnames(indel) <- indel_types
      }

      if('indel' %in% lsqnonneg.types){
         stop('Indel least-squares fit signatures have not been implemented in this version')
      }

      if('indel' %in% rel.types){
         indel <- indel/rowSums(indel)
      }

      out[['indel']] <- indel
   }

   if(!is.null(sv)){
      if('sv' %in% simplify.types){
         sv_types <- c('DEL', 'DUP', 'INV', 'TRA')
         sv_split <- splitDfRegex(sv_types, df = sv)
         sv <- do.call(cbind, lapply(sv_split, rowSums))
         colnames(sv) <- paste0('SV.',sv_types)
      }

      if('sv' %in% lsqnonneg.types){
         sv <- t(apply(contexts$sv, 1, function(i){
            fitToSignatures(SV_SIGNATURE_PROFILES, i)$x
         }))
         colnames(sv) <- paste0('SV',1:6)
      }

      if('sv' %in% rel.types){
         sv <- sv/rowSums(sv)
      }

      out[['sv']] <- sv
   }

   do.call(cbind, unname(out))
}


