#' Read vcf into R as data frame
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

   ## Scan for the header line
   con  <- file(vcf.file, open = "r")
   line_num <- 0
   while(length(line <- readLines(con, n=1, warn=F)) > 0) {
      line_num <- line_num + 1
      #if(grepl('^#CHROM',line)){ print(line) }
      #if(line_num==100){ break }
      if(!grepl('^##',line)){ break }
   }
   close(con)

   vcf <- read.delim(
      vcf.file, skip=line_num-1,
      check.names=F, stringsAsFactors=F,
      colClasses='character'
   )
   vcf$POS <- as.integer(vcf$POS)

   ## Remove '#' from header line
   colnames(vcf) <- sub('^#','',colnames(vcf))

   ## Select desired columns
   if(!is.null(fields)){
      vcf <- vcf[,fields,drop=F]
   }

   return(vcf)
}



####################################################################################################
#' Get values from INFO field
#'
#' @param v A character vector of the INFO field, with each item being a line of the INFO field
#' @param keys A character vector of the names of the INFO field values to retrieve
#'
#' @return A character matrix containing the key names and corresponding values
#' @export
getInfoValues <- function(v, keys){
   #v=vcf$info
   l <- strsplit(v,';')
   l <- lapply(l, function(i){
      do.call(rbind,strsplit(i,'='))
   })

   out <- do.call(rbind,lapply(l, function(i){
      #i=l[[1]]
      v <- i[match(keys, i[,1]),2]
      if(length(v)==0){
         return(rep(NA,length(keys)))
      }
      return(v)
   }))

   colnames(out) <- keys

   return(as.data.frame(out, stringsAsFactors=F))
}
