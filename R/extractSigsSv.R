#' Extract structural variant signatures
#'
#' @param vcf.file Path to the vcf file
#' @param output Output the absolute signature contributions (default, 'signatures'), or the SV type/length contexts
#' ('contexts')
#' @param sample.name If a character is provided, the header for the output matrix will be named to this. If none is
#' provided, the basename of the vcf file will be used.
#' @param sv.len.cutoffs SV length cutoff intervals as a numeric vector.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are the mutation
#' contexts and the columns are  the mutational signatures.
#'
#' @return A 1-column matrix containing: (if output = 'signatures') the absolute signature contributions (i.e. the number of mutations contributing to
#' each mutational signature) , or (if output = 'contexts') the mutation contexts.
#' @export

extractSigsSv <- function(vcf.file, output = 'signatures', sample.name = NULL, sv.caller = 'manta',
                          sv.len.cutoffs = c(10^3, 10^4, 10^5, 10^6, 10^7, Inf),
                          signature.profiles = SV_SIGNATURE_PROFILES, ...){

   variants <- variantsFromVcf(vcf.file, mode = 'sv', sv.caller = sv.caller, ...)

   ## Initialize table of SV type/length bins
   sv_types <- c('DEL','DUP','INV') ## INS ignored. TRA/BND dealt with in a later step

   sv_contexts <- data.frame(
      sv_type = rep( sv_types, each = length(sv.len.cutoffs)-1 ),
      lower_cutoff = rep( sv.len.cutoffs[-length(sv.len.cutoffs)], length(sv_types) ),
      upper_cutoff = rep( sv.len.cutoffs[-1], length(sv_types) ),

      stringsAsFactors = F
   )

   ## Deal with empty vcfs
   if(!is.data.frame(variants) && is.na(variants)){
      context_counts <- rep(0, nrow(sv_contexts)+1)
   }

   else {
      ## Count context occurrences
      context_counts <- unlist(lapply(1:nrow(sv_contexts), function(i){
         row <- as.list(sv_contexts[i,])

         variants_ss <- variants[
            variants$sv_type == row$sv_type
            & variants$sv_len >= row$lower_cutoff
            & variants$sv_len < row$upper_cutoff
            ,]

         nrow(variants_ss)
      }))

      ## Count context occurrences for translocations
      translocation_counts <- nrow(variants[variants$sv_type == 'BND' | variants$sv_type == 'TRA',])

      if(sv.caller == 'manta'){ ## manta reports translocations twice (origin/destination)
         translocation_counts <- translocation_counts/2
      }

      context_counts <- c(context_counts,translocation_counts)
   }

   ## Create context names
   names(context_counts) <- c(
      str_replace_all(paste(
            sv_contexts$sv_type,
            formatC(sv_contexts$lower_cutoff, format = 'e', digits = 0),
            formatC(sv_contexts$upper_cutoff, format = 'e', digits = 0),
            'bp', sep = '_'
         ),'[+]',''),

      'TRA'
   )

   if(output == 'contexts'){
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){

      if(nrow(signature.profiles) != length(context_counts)){
         stop('The number of contexts in the signature profile matrix != the number of contexts in the context count vector.\n
              Check that the provided cutoffs sv.len.cutoffs also exists in signature.profiles')
      }

      ## Least squares fitting
      out <- fitToSignatures(signature.profiles, context_counts)$x
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }

   if(is.null(sample.name)){ colnames(out) <- basename(vcf.file) }
   else { colnames(out) <- sample.name }

   return(out)
}

# #========= Testing =========#
# XXXXXXXX <- list(
#    gridss = extractSigsSv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mutSigExtractor/R_source/gridss_output/XXXXXXXX.purple.sv.vcf.gz',
#                           sv.caller = 'gridss', output = 'contexts'),
#    bpi = extractSigsSv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170812_HMFregXXXXXXXX/XXXXXXXX.vcf.gz',
#                        sv.caller = 'manta', output = 'contexts')
# )
#
# XXXXXXXX <- list(
#    gridss = extractSigsSv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mutSigExtractor/R_source/gridss_output/XXXXXXXX.purple.sv.vcf.gz',
#                           sv.caller = 'gridss', output = 'signatures'),
#    bpi = extractSigsSv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170812_HMFregXXXXXXXX/XXXXXXXX.vcf.gz',
#                        sv.caller = 'manta', output = 'signatures')
# )
#
# do.call(cbind, XXXXXXXX)
# do.call(cbind, XXXXXXXX)
#
# extractSigsIndel('/Users/lnguyen//hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data/170812_HMFregXXXXXXXX/XXXXXXXX.vcf.gz')
