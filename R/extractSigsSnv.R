#' Extract single nucleotide variant signatures
#'
#' @description Will output a 1-column matrix containing: (if output = 'signatures') the absolute
#' signature contributions (i.e. the number of mutations contributing to each mutational signature),
#' or (if output = 'contexts') the mutation contexts.
#'
#' To elaborate the signatures used are the 30 COSMIC signatures (https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt),
#' which are derived from the 96-trinucleotide mutation contexts.
#'
#' @param vcf.file Path to the vcf file
#' @param output Output the absolute signature contributions (default, 'signatures'), or the
#' 96-trinucleotide contexts ('contexts')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param ref.genome A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#'
#' @return A 1-column matrix
#' @export
extractSigsSnv <- function(vcf.file, output = 'signatures', sample.name = NULL,
                           ref.genome = DEFAULT_GENOME, signature.profiles = SNV_SIGNATURE_PROFILES, ...){

   #vcf.file = '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/R_source/call_snv_signatures/EGAZ00001222791_EGAZ00001222814_post_processed.vcf.gz'
   #output = 'contexts'
   variants <- variantsFromVcf(vcf.file, mode = 'snv', ref.genome, ...)
   #variants <- variantsFromVcf(vcf.file, mode = 'snv', ref.genome)

   ## Deal with empty vcfs
   if(!is.data.frame(variants) && is.na(variants)){
      subs_context_counts <- rep(0, length(SUBS_CONTEXTS_96))
      names(subs_context_counts) <- SUBS_CONTEXTS_96
   }

   else {
      ## Get trinucleotide contexts that don't conform to C>N or T>N
      select_opp_types <- which(!(variants$substitution %in% SUBSTITUTIONS))

      ## Reverse complement for non-conforming contexts
      variants[select_opp_types,'tri_context'] <- reverse(
         chartr('ATGC', 'TACG', variants[select_opp_types,'tri_context'])
      )

      ## Single nt can simply be complemented (no need for reverse)
      variants[select_opp_types,'substitution'] <- chartr('ATGC', 'TACG', variants[select_opp_types,'substitution'])

      ## Convert trinucleotide context to substitution context
      variants$subs_context <- paste0(
         substr(variants$tri_context, 1, 1),
         '[', variants$substitution, ']',
         substr(variants$tri_context, 3, 3)
      )

      ## Count context occurrences. Fill unfound contexts with 0
      subs_context_counts <- table(variants$subs_context)

      subs_context_counts <- unlist(lapply(SUBS_CONTEXTS_96, function(i){
         counts <- subs_context_counts[i]
         if(is.na(counts)){ counts <- 0 }
         return(counts)
      }))

      names(subs_context_counts) <- SUBS_CONTEXTS_96
   }

   ## Export
   if(output == 'contexts'){
      out <- as.matrix(subs_context_counts)

   } else if(output == 'signatures'){
      ## Least squares fitting
      out <- fitToSignatures(signature.profiles, subs_context_counts)$x
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }

   if(is.null(sample.name)){ colnames(out) <- basename(vcf.file) }
   else { colnames(out) <- sample.name }

   return(out)
}
