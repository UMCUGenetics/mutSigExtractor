#' Extract trinucleotide contexts
#'
#' @param bed A dataframe containing the columns: chrom, pos, ref, alt
#' @param ref.genome A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param verbose Print progress messages?
#'
#' @return A dataframe in the same structure as a bed file
#' @export
getContextsSnv <- function(bed, ref.genome=DEFAULT_GENOME, verbose=F){

   bed_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(bed)[1:4], bed_colnames))){
      warning("colnames(bed)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(bed)[1:4] <- bed_colnames
   }

   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   bed <- bed[!grepl(',',bed$alt),]

   if(verbose){ message('Subsetting for SNVs...') }
   bed <- bed[nchar(bed$ref)==1 & nchar(bed$alt)==1,]
   if(nrow(bed)==0){
      warning('No variants remained after subsetting for SNVs. Returning NA')
      return(NA)
   }

   if(verbose){ message('Converting chrom name style to style in ref.genome...') }
   seqlevelsStyle(bed$chrom) <- seqlevelsStyle(eval(parse(text=ref.genome)))

   if(verbose){ message('Returning SNV trinucleotide contexts...') }
   out <- data.frame(
      substitution = paste0(bed$ref,'>',bed$alt),
      tri_context = getSeq(
         x = eval(parse(text=ref.genome)),
         names = bed$chrom,
         start = bed$pos - 1,
         end = bed$pos + 1,
         as.character=T
      ),
      stringsAsFactors=F
   )

   return(out)
}

####################################################################################################

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
#' @param bed A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to vcf.file
#' @param output Output the absolute signature contributions (default, 'signatures'), or the
#' 96-trinucleotide contexts ('contexts')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param ref.genome A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param verbose Print progress messages?
#'
#' @return A 1-column matrix
#' @export
extractSigsSnv <- function(
   vcf.file=NULL, bed=NULL, output='signatures', sample.name=NULL,
   ref.genome=DEFAULT_GENOME, signature.profiles=SNV_SIGNATURE_PROFILES,
   verbose=F, ...
){

   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      bed <- variantsFromVcf(vcf.file, mode='snv_indel', ref.genome=ref.genome, verbose=verbose, ...)
   }
   df <- getContextsSnv(bed, ref.genome=ref.genome, verbose=verbose)

   if(verbose){ message('Initializing SNV signature output vector...') }
   context_counts <- structure(
      rep(0, length(SUBS_CONTEXTS_96)),
      names=SUBS_CONTEXTS_96
   )

   if(is.data.frame(df)){ ## Don't process empty vcfs (df==NA if empty)
      if(verbose){ message('Converting trinucleotide contexts to substitution contexts...') }
      ## Get trinucleotide contexts that don't conform to C>N or T>N
      select_opp_types <- which(!(df$substitution %in% SUBSTITUTIONS))

      ## Reverse complement for non-conforming contexts
      df[select_opp_types,'tri_context'] <- reverse(
         chartr('ATGC', 'TACG', df[select_opp_types,'tri_context'])
      )

      ## Single nt can simply be complemented (no need to reverse)
      df[select_opp_types,'substitution'] <- chartr('ATGC', 'TACG', df[select_opp_types,'substitution'])

      ## Convert trinucleotide context to substitution context
      df$subs_context <- paste0(
         substr(df$tri_context, 1, 1),
         '[', df$substitution, ']',
         substr(df$tri_context, 3, 3)
      )

      ## Count context occurrences. Fill found contexts into context_counts
      if(verbose){ message('Counting substitution context occurrences...') }
      context_counts_new <- table(df$subs_context)
      context_counts[names(context_counts_new)] <- context_counts_new
   }

   ## Export
   if(output == 'contexts'){
      if(verbose){ message('Returning context counts...') }
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){
      if(verbose){ message('Returning absolute signature contributions...') }
      ## Least squares fitting
      out <- fitToSignatures(signature.profiles, context_counts)$x
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }
   if(verbose & !is.data.frame(df)){ warning("Input to extractSigsSnv() contained no variants. Returning dataframe of 0's") }

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}
