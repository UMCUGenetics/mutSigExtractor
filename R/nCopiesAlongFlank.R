#' Calculate the number of copies of the indel sequence are present in the flank sequence
#'
#' @description Helper function for extractSigsIndel(). Scans along the flanking sequence for copies of the indel sequence
#' using a sliding window with a length equal to the length of the indel sequence. The sliding window jumps in increments equal
#' to the indel length. Scanning stops prematurely when the indel sequence does not match the one in the sliding window.
#' The indel sequence itself counts as one copy.
#'
#' @param indel.seq The indel sequence as a character
#' @param flank.seq The flanking sequence as a character
#'
#' @return An integer stating the number of copies found in the flanking sequence
#' @export
nCopiesAlongFlank <- function(indel.seq, flank.seq){

   # indel.seq = "T"
   # flank.seq = "TTTTGCG"

   indel_length <- nchar(indel.seq)
   rail_seq <- paste0(indel.seq, flank.seq)

   count <- 0
   seq_window <- substr(rail_seq, 1, indel_length)

   while(indel.seq == seq_window){
      count <- count + 1
      seq_window <- substr(rail_seq,
                           count*indel_length + 1,
                           count*indel_length + indel_length)
   }

   return(count)
}
