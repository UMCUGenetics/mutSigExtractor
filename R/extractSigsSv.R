#' Extract structural variant signatures
#'
#' @description Will return a 1-column matrix containing: (if output = 'signatures') the absolute
#' signature contributions (i.e. the number of mutations contributing to each mutational signature),
#' or (if output = 'contexts') the mutation contexts.
#'
#' To elaborate, the 6 SV signatures used are those described in this paper: https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip,
#' in Supplementary.Table.21.Signatures.v3.xlsx. These are derived from mutation contexts composed
#' of SV type/length.
#'
#' Note that the probabilities of the clustered and non-clustered rearrangements in the signature
#' profile have been combined. In other words, whether the rearrangements were
#' clustered/non-clustered were not considered.
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe with the columns: sv_type, sv_len. sv_type can be DEL, DUP, INV, TRA, BND.
#' Note that for TRA and BND, sv_len will be ignored. Alternative input option to vcf.file
#' @param output Output the absolute signature contributions (default, 'signatures'), or the SV
#' type/length contexts ('contexts')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param half.tra.counts Divide translocation counts by 2?
#' @param sv.caller Can be 'manta' or 'gridss'
#' @param sv.len.cutoffs SV length cutoff intervals as a numeric vector.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param verbose Print progress messages?
#'
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export

extractSigsSv <- function(
   vcf.file=NULL, df=NULL, output='signatures', sample.name=NULL,
   sv.caller='gridss', half.tra.counts=F,
   sv.len.cutoffs=if(output=='signatures'){ c(10^c(3:7), Inf) } else { c(0, 10^c(3:7), Inf) },
   signature.profiles=SV_SIGNATURE_PROFILES,
   verbose=F, ...
){
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, mode='sv', sv.caller=sv.caller, verbose=verbose, ...)
   } else if(is.data.frame(df)){
      colnames(df) <- c('sv_type','sv_len')
      half.tra.counts <- F ## If providing dataframe as input default to 'manta'.
   } else if(!is.na(df)) {
      stop('Please specify either vcf.file or df as input')
   }

   if(verbose){ message('Creating SV type/length lookup table...') }
   sv_types <- c('DEL','DUP','INV') ## INS ignored. TRA/BND dealt with in a later step

   sv_contexts <- data.frame(
      sv_type = rep( sv_types, each = length(sv.len.cutoffs)-1 ),
      lower_cutoff = rep( sv.len.cutoffs[-length(sv.len.cutoffs)], length(sv_types) ),
      upper_cutoff = rep( sv.len.cutoffs[-1], length(sv_types) ),

      stringsAsFactors = F
   )

   sv_contexts$name <- with(sv_contexts,{
      v <- paste(
         sv_type,
         formatC(lower_cutoff, format = 'e', digits = 0),
         formatC(upper_cutoff, format = 'e', digits = 0),
         'bp',sep='_'
      )
      gsub('[+]','',v)
   })

   ## Deal with empty vcfs
   if(!is.data.frame(df) && is.na(df)){
      context_counts <- rep(0, nrow(sv_contexts)+1)
   }

   else {
      if(verbose){ message('Counting DEL, DUP, and INV context occurrences...') }
      context_counts <- unlist(lapply(1:nrow(sv_contexts), function(i){
         row <- sv_contexts[i,]
         variants_ss <- df[
            df$sv_type == row$sv_type
            & df$sv_len >= row$lower_cutoff
            & df$sv_len < row$upper_cutoff
            ,]

         nrow(variants_ss)
      }))

      ## Count context occurrences for translocations
      if(verbose){ message('Counting TRA occurrences...') }
      translocation_counts <- nrow(df[df$sv_type == 'BND' | df$sv_type == 'TRA',])

      if(sv.caller=='manta' & half.tra.counts){ ## manta reports translocations twice (origin/destination)
         if(verbose){ message('Halving TRA counts...') }
         translocation_counts <- translocation_counts/2
      }

      context_counts <- c(context_counts,translocation_counts)
   }

   ## Assign context names
   names(context_counts) <- c(sv_contexts$name,'TRA')

   if(output == 'contexts'){
      if(verbose){ message('Returning SV contexts...') }
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){
      if(verbose){ message('Returning SV signatures...') }

      if(nrow(signature.profiles) != length(context_counts)){
         stop('The number of contexts in the signature profile matrix != the number of contexts in the context count vector.\n
              Check that the provided cutoffs sv.len.cutoffs also exists in signature.profiles')
      }

      ## Least squares fitting
      out <- fitToSignatures(signature.profiles, context_counts)$x
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}
