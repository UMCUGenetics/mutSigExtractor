#' Get the start/end positions for the left/right flanks of an indel
#'
#' @description Helper function for extractSigsIndel(), getting flank start/end positions first, then retrieving the sequences
#' using getSeq with a vectorized input improves speed significantly.
#'
#' @param chrom A character stating the chromosome
#' @param pos The REF position as an integer
#' @param indel.len The length of the indel as an integer
#' @param indel.type A character stating whether the indel is an insertion ('ins') or deletion ('del')
#' @param n.indel.lengths.l The length of the flanking sequence to return (measured in indel lengths)
#' @param n.indel.lengths.r See n.indel.lengths.l
#'
#' @return A vector of the left flank start/end position and right flank start/end position
#' @export
indelSeqFlanksStartEnd <- function(chrom, pos, indel.len, indel.type,
                                   n.indel.lengths.l = 1, n.indel.lengths.r = 1){

   if(indel.type == 'del'){
      r_flank <- c(
         r_start = pos + 1 + indel.len,
         r_end = pos + indel.len + indel.len*n.indel.lengths.r
      )

   } else if(indel.type == 'ins'){
      r_flank <- c(
         r_start = pos + 1,
         r_end = pos + indel.len*n.indel.lengths.r
      )
   }

   l_flank <- c(
      l_start = pos - indel.len*n.indel.lengths.l + 1,
      l_end = pos
   )

   out <- l_flank
   out <- c(out, r_flank)

   return(out)
}
