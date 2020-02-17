#' Extract doublet base substitution contexts
#'
#' @param df A dataframe containing the columns: chrom, pos, ref, alt
#' @param ref.genome A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param verbose Print progress messages?
#'
#' @return A dataframe in the same structure as a bed file with an extra column stating the context
#' of each variant
#' @export
getContextsDbs <- function(df, ref.genome=DEFAULT_GENOME, verbose=F){

   df_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(df)[1:4], df_colnames))){
      warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(df)[1:4] <- df_colnames
   }

   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]

   if(verbose){ message('Converting chrom name style to style in ref.genome...') }
   seqlevelsStyle(df$chrom) <- seqlevelsStyle(eval(parse(text=ref.genome)))

   if(verbose){ message('Subsetting for DBSs...') }
   df <- df[nchar(df$ref)==2 & nchar(df$alt)==2,]
   if(nrow(df)==0){
      warning('No variants remained after subsetting for DBSs. Returning NA')
      return(NA)
   }

   ## Check for weird nucleotides
   which_weird_nt <- sort(unique(c(
      grep('[^ACTG]',df$ref),
      grep('[^ACTG]',df$alt)
   )))

   if(length(which_weird_nt)>0){
      warning(
         length(which_weird_nt),
         ' variants containing nucleotides other than A,T,C,G were removed (rows: ',
         paste(which_weird_nt, collapse=', '), ')'
      )
      df <- df[-which_weird_nt,]
      if(nrow(df)==0){
         warning('No variants remained after removing weird nucleotides. Returning NA')
         return(NA)
      }
   }

   df$context <- paste0(df$ref,'>',df$alt)

   ## Get reverse complement (from lookup table) where applicable
   df$index <- 1:nrow(df) ## Create index to maintain original row order
   df_split <- split(df, df$context %in% DBS_TYPES$context)

   df_split[['FALSE']] <- within(df_split[['FALSE']],{
      context <- DBS_TYPES$context[ match(context, DBS_TYPES$context_rev_comp) ]
   })

   df <- rbind(df_split[['FALSE']],df_split[['TRUE']]); rm(df_split)
   df <- df[order(df$index),]; df$index <- NULL ## Restore original row order

   return(df)
}


####################################################################################################

#' Extract doublet substitution signatures
#'
#' @description Will output a 1-column matrix containing: (if output = 'signatures') the absolute
#' signature contributions (i.e. the number of mutations contributing to each mutational signature),
#' or (if output = 'contexts') the mutation contexts.
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to vcf.file
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
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export
extractSigsDbs <- function(
   vcf.file=NULL, df=NULL, output='contexts', sample.name=NULL,
   ref.genome=DEFAULT_GENOME, signature.profiles=DBS_SIGNATURE_PROFILES,
   verbose=F, ...
){

   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, mode='snv_indel', ref.genome=ref.genome, verbose=verbose, ...)
   }
   df <- getContextsDbs(df, ref.genome=ref.genome, verbose=verbose)

   ## Check for weird nucleotides
   which_weird_nt <- sort(unique(c(
      grep('[^ACTG>]',df$substitution),
      grep('[^ACTG]',df$tri_context)
   )))

   if(length(which_weird_nt)>0){
      warning(
         length(which_weird_nt),
         ' variants containing nucleotides other than A,T,C,G were removed (rows: ',
         paste(which_weird_nt, collapse=', '), ')'
      )
      df <- df[-which_weird_nt,]
   }

   if(verbose){ message('Counting DBS context occurrences...') }
   context_counts <- structure(
      rep(0, nrow(DBS_TYPES)),
      names=DBS_TYPES$context
   )

   if(is.data.frame(df)){
      tab <- table(df$context)
      context_counts[names(tab)] <- tab
   }


   if(output == 'contexts'){
      if(verbose){ message('Returning DBS context counts...') }
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){
      if(verbose){ message('Returning DBS signature contributions...') }
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


