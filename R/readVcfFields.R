#' A basic vcf reader
#'
#' @param vcf.file Path to the vcf file
#' @param fields A character or integer vector indicating which columns to keep.
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' readVcfFields('/path/to/vcf', fields=c('CHROM','POS','REF','ALT'))
#' readVcfFields('/path/to/vcf', fields=c(1,2,4,5))

readVcfFields <- function(vcf.file, fields=NULL){

   ## Remove all header lines, then read the vcf
   clean_lines <- sub('##.*','', readLines(vcf.file))
   vcf <- read.delim(text=paste(clean_lines, collapse = "\n"), check.names=F, stringsAsFactors=F)

   ## Remove '#' from header line
   colnames(vcf) <- sub('^#','',colnames(vcf))

   ## Select desired columns
   if(!is.null(fields)){
      vcf <- vcf[,fields,drop=F]
   }

   if(nrow(vcf)==0){
      warning('VCF contains no rows. Returning NA')
      stop(return(NA))
   }

   return(vcf)
}

