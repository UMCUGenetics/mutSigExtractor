#' Subset for SNVs, DBSs, or indels
#'
#' @description Standard variant filtering and subsetting for variant type used to clean the input
#' dataframe for extractSigs{Snv,Indel,Dbs}() functions
#'
#' @param df A dataframe containing the columns chrom, pos, ref, alt
#' @param type Can be 'snv', 'indel', or 'dbs'
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param verbose Show progress messages?
#'
#' @return The filtered dataframe, or an empty dataframe if no variants exist/remain
#' @export
#'
subsetSmnvs <- function(df, type, ref.genome=DEFAULT_GENOME, verbose=F){

   if(length(type)!=1 | !(type %in% c('snv','indel','dbs'))){
      stop("`type` must be one of the following values: 'snv','indel','dbs'")
   }

   ## Filter variants --------------------------------
   if(nrow(df)==0){
      return(data.frame())
   }

   df_colnames <- c('chrom','pos','ref','alt')
   if(!(identical(colnames(df)[1:4], df_colnames))){
      warning("colnames(df)[1:4] != c('chrom','pos','ref','alt'). Assuming first 4 columns are these columns")
      colnames(df)[1:4] <- df_colnames
   }

   if(verbose){ message('Removing rows with multiple ALT sequences...') }
   df <- df[!grepl(',',df$alt),]

   if(nrow(df)==0){
      warning('No variants remained after subsetting for variants with one ALT sequence. Returning empty dataframe')
      return(data.frame())
   }

   #if(verbose){ message('Converting chrom name style to style in ref.genome...') }
   GenomeInfoDb::seqlevelsStyle(df$chrom)<- GenomeInfoDb::seqlevelsStyle(ref.genome)

   ## Select variant type --------------------------------
   df$ref_len <- nchar(df$ref)
   df$alt_len <- nchar(df$alt)

   if(type=='snv'){
      if(verbose){ message('Subsetting for SNVs...') }
      df <- df[df$ref_len==1 & df$alt_len==1,]
      if(nrow(df)==0){
         warning('No variants remained after subsetting for SNVs. Returning empty dataframe')
         return(data.frame())
      }
   }

   if(type=='dbs'){
      if(verbose){ message('Subsetting for DBSs...') }
      df <- df[df$ref_len==2 & df$alt_len==2,]
      if(nrow(df)==0){
         warning('No variants remained after subsetting for DBSs. Returning empty dataframe')
         return(data.frame())
      }
   }

   if(type=='indel'){
      if(verbose){ message('Determining indel type...') }
      ## Calc sequence lengths


      ## Remove snvs
      df <- df[!(df$ref_len==1 & df$alt_len==1),]
      if(nrow(df)==0){
         warning('No variants remained after subsetting for indels. Returning empty data.frame')
         return(data.frame())
      }
   }

   return(df)
}
