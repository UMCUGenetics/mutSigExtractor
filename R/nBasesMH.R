#' Calculate the number of bases in microhomology
#'
#' @description Helper function for extractSigsIndel(). Scans a maximum of 1 indel length from 5' to 3' in the flanking sequence to
#' detect identical bases. Stops prematurely if a non-identical base is found.
#'
#' @param indel.seq The indel sequence as a character
#' @param flank.seq The flanking sequence as a character
#'
#' @return An integer stating the number of bases in microhomology
#' @export
nBasesMH <- function(indel.seq, flank.seq, indel.len){

   #indel.sequence = "CTA"
   #flank.sequence = "C"

   indel_len <- nchar(indel.seq)

   indel.seq <- unlist(str_split(indel.seq, ''))
   flank.seq <- unlist(str_split(flank.seq, ''))[1:indel_len]

   n_bases <- 0
   for(i in 1:length(indel.seq)){
      if(indel.seq[i] != flank.seq[i]){
         break
      } else {
         n_bases <- n_bases + 1
      }
   }
   return(n_bases)
}
