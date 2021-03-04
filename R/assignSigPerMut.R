#' Assign a mutational signature to each mutation in a vcf or bed-like dataframe
#'
#' @description This function performs the following steps to assign a signature to each mutation
#' in a vcf or bed-like dataframe:
#' (1) Extract contexts and count occurrences
#' (2) Calculate signature contributions by fitting contexts counts to reference signature profiles
#' (3) Multiply each reference signature profile by the signature contribution vector. This gives
#' the probability of each context to be assigned to a signature
#' (4) Assign each context a signature based on the maximum probability signature
#' (5) Assign each mutation a signature based on its context
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
#' vcf.file
#' @param mode Can be 'snv','sbs','indel' or 'dbs'. `mode` defines which `extractSigs*()` function
#' to use. Additionally, if `signature.profiles` is unspecified, the correct signature profile
#' matrix will be automatically chosen
#' @param output Can be 'df' (a dataframe with chrom, pos, ref, alt, context, assigned_sig, sig_prob),
#' 'df.compact' (a dataframe with context, assigned_sig, sig_prob), or 'vector' (a vector whose
#' names are assigned_sig, and values are sig_prob).
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#' @param args.extract.sigs Args that can be passed to the `extractSigs*()` functions
#' @param fit.method Can be 'lsq' or 'strict'. Method for fitting context counts to signature
#' profiles. See documentation for `fitToSignatures()` for more details
#' @param args.fit Args that can be passed to the `fitToSignatures*()` functions
#' @param verbose Print progress messages?
#'
#' @return See `output`
#' @export
#'
assignSigPerMut <- function(
   vcf.file=NULL, df=NULL, mode='snv', output='df',
   signature.profiles=NULL,
   args.extract.sigs=list(vcf.filter='PASS'),
   fit.method='strict', args.fit=NULL,
   verbose=T
){

   if(F){
      vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/200706_HMFregXXXXXXXX/purple/XXXXXXXX.purple.somatic.vcf.gz'
      mode='snv'
      output='df'
      signature.profiles=SBS_SIGNATURE_PROFILES_V3
      fit.method='strict'
      args.fit=NULL
      args.extract.sigs=list(vcf.filter='PASS')
      verbose=T
   }

   ## Checks --------------------------------
   if(!(mode %in% c('snv','sbs','indel','dbs'))){
      stop("`mode` must be one of the following: 'snv','sbs','indel', 'dbs'")
   }

   if(!(fit.method %in% c('lsq','strict'))){
      stop("`fit.method` must be one of the following: 'lsq','strict'")
   }

   if(!(output %in% c('df','df.compact','vector'))){
      stop("`output` must be one of the following: 'df','df.compact','vector'")
   }

   ## --------------------------------
   if(verbose){ message('Extracting ',mode,' contexts...') }
   f_extract_sigs <- switch(
      mode,
      snv=extractSigsSnv,
      sbs=extractSigsSnv,
      dbs=extractSigsDbs,
      indel=extractSigsIndel
   )

   ## Oligatory args
   args_extract_sigs <- list(vcf.file=NULL, df=NULL, output='df')
   if(mode=='indel'){ args_extract_sigs <- c(args_extract_sigs, method='PCAWG') }

   ## Optional args
   args_extract_sigs <- c(args_extract_sigs, args.extract.sigs)

   ## Rm duplicate args
   args_extract_sigs <- args_extract_sigs[ !duplicated(names(args_extract_sigs)) ]

   if(!is.null(vcf.file)){
      args_extract_sigs$vcf.file <- vcf.file
   } else if(!is.null(df)){
      args_extract_sigs$df <- df
   } else {
      stop('Input must be specified to `vcf.file` or `df`')
   }

   df <- do.call(f_extract_sigs, args_extract_sigs)
   rm(f_extract_sigs, args_extract_sigs)

   ## Select only relevant columns
   df <- df[,c('chrom','pos','ref','alt','context')]

   ## --------------------------------
   if(verbose){ message('Counting mutation contexts...') }
   context_counts_pre <- table(df$context)
   context_counts <- structure(
      rep(0, length(levels(df$context))),
      names=levels(df$context)
   )
   context_counts[names(context_counts_pre)] <- context_counts_pre

   ## --------------------------------
   if(verbose){ message('Fitting context counts to signature profiles...') }
   if(!is.null(signature.profiles)){
      sig_profiles <- signature.profiles
   } else {
      sig_profiles <- switch(
         mode,
         snv=SBS_SIGNATURE_PROFILES_V3,
         sbs=SBS_SIGNATURE_PROFILES_V3,
         dbs=DBS_SIGNATURE_PROFILES,
         indel=INDEL_SIGNATURE_PROFILES
      )
   }
   #sig_profiles <- sig_profiles / colSums(sig_profiles) ## Force total prob to be 1

   ## Force context names and order to be the same in sig profile matrix
   if(!all(context_counts %in% rownames(sig_profiles))){
      stop('\nExtracted contexts do not match `rownames(signature.profiles)`\nMaybe the wrong `mode` was selected?')
   }
   sig_profiles <- sig_profiles[names(context_counts),]

   f_fit <- switch(
      fit.method,
      lsq=fitToSignaturesFast,
      strict=fitToSignaturesFastStrict
   )

   args_fit <- c(
      list(mut.context.counts=context_counts, signature.profiles=sig_profiles),
      args.fit
   )
   args_fit <- args_fit[ !duplicated(names(args_fit)) ]

   sig_contrib <- do.call(f_fit, args_fit)
   rm(f_fit, args_fit)

   ## --------------------------------
   if(verbose){ message('Calculating signature probabilities per mutation...') }

   ## Adjust signature context probabilities based on signature contributions in sample
   sig_profiles_sample <- t(apply(sig_profiles, 1, function(i){ i * sig_contrib }))
   colnames(sig_profiles_sample) <- colnames(sig_profiles)
   sig_profiles_sample <- sig_profiles_sample / rowSums(sig_profiles_sample)

   ## Get max probability signature per context
   context_sig_assignment <- data.frame(
      assigned_sig = colnames(sig_profiles_sample)[ max.col(sig_profiles_sample) ],
      sig_prob = apply(sig_profiles_sample,1,max),
      row.names=rownames(sig_profiles_sample)
   )

   ## Assign each mutation a signature based on its context
   df <- data.frame(
      df,
      context_sig_assignment[df$context,]
   )

   ## --------------------------------
   if(verbose){ message('Returning output...') }

   if(output=='df'){ return(df) }

   if(output=='df.compact') {
      return( df[,c('context','assigned_sig','sig_prob')] )
   }

   structure(df$sig_prob, names=df$assigned_sig)
}
