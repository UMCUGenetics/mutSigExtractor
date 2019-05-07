#' Extract indel signatures
#'
#' @param description Will return a 1-column matrix containing the absolute indel signature
#' contributions (i.e. the number of mutations contributing to each mutational signature). The
#' signatures used are insertions/deletions within repeat regions (ins.rep, del.rep),
#' insertions/deletions with flanking microhomology (ins.mh, del.mh), and insertions/deletions
#' which don't fall under the previous 2 categories (ins.none, del.none). Each category is further
#' stratified by the length of the indel.
#'
#' @param vcf.file Path to the vcf file
#' @param sample.name If a character is provided, the header for the output matrix will be named to this. If none is
#' provided, the basename of the vcf file will be used.
#' @param ref.genome A character naming the BSgenome reference genome. Default is "BSgenome.Hsapiens.UCSC.hg19". If another
#' reference genome is indicated, it will also need to be installed.
#' @param indel.len.cap Specifies the max indel sequence length to consider when counting 'repeat' and 'none' contexts.
#' Counts of longer indels will simply be binned to the counts of contexts at the max indel sequence length.
#' @param n.bases.mh.cap Specifies the max bases in microhomology to consider when counting repeat and microhomology
#' contexts. Counts of longer indels will simply be binned to the counts of contexts at the max indel sequence length.
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to variantsFromVcf()
#'
#' @return A 1-column matrix
#' @export

extractSigsIndel <- function(
   vcf.file, sample.name=NULL, ref.genome=DEFAULT_GENOME,
   indel.len.cap=5, n.bases.mh.cap=5, verbose=F, ...
){
   variants <- variantsFromVcf(vcf.file, mode = 'indel', ref.genome, verbose=verbose, ...)
   # variants <- variantsFromVcf(
   #    vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/ICGC/vcf/BRCA-EU/snv_indel/PD10010_snv_indel.vcf.gz',
   #    mode = 'indel', ref.genome
   # )

   # variants <- variantsFromVcf(
   #    '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/ICGC/vcf/OV-AU/snv_indel/AOCS-001_snv_indel.vcf.gz'
   #    ,mode='indel',get.other.indel.allele=T
   # )

   if(verbose){ message('Initializing indel signature output vector...') }
   indel_sig_names <- c(
      paste0('del.rep.len.', 1:indel.len.cap),
      paste0('ins.rep.len.', 1:indel.len.cap),
      paste0('del.mh.bimh.', 1:n.bases.mh.cap),
      paste0('ins.mh.bimh.', 1:n.bases.mh.cap),
      paste0('del.none.len.', 1:indel.len.cap),
      paste0('ins.none.len.', 1:indel.len.cap)
   )
   indel_sigs <- structure(rep(0,length(indel_sig_names)), names=indel_sig_names)

   if(is.data.frame(variants) && !is.na(variants)){ ## Deal with empty vcfs

      #--------- Pre-calculations for repeat and microhomology contexts ---------#
      if(verbose){ message('Determining the start/end positions for the left/right flanks of each indel...') }
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
      indelSeqFlanksStartEnd <- function(chrom, pos, indel.len, indel.type, n.indel.lengths.l=1, n.indel.lengths.r=1){
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

      ## Getting flank start/end positions first, then retrieving the sequences using getSeq with a
      ## vectorized input improves speed significantly.
      ## Cap n.indel.lengths.r to 3 (length of 3' (right-hand side) sequence to retrieve),
      ## since the max n_copies_along_flank condition used below caps at >=2.
      ## This improves speed significantly
      flanks_start_end <- with(variants,{
         do.call(rbind, Map(
            f = indelSeqFlanksStartEnd,
            chrom, pos, indel_len, indel_type,
            n.indel.lengths.r=3
         ))
      })

      if(verbose){ message('Retrieving flanking sequences...') }
      l_flank <- getSeq(
         x = eval(parse(text=ref.genome)),
         names = variants$chrom,
         start = flanks_start_end[,'l_start'],
         end = flanks_start_end[,'l_end'],
         as.character = T
      )

      r_flank <- getSeq(
         x = eval(parse(text=ref.genome)),
         names = variants$chrom,
         start = flanks_start_end[,'r_start'],
         end = flanks_start_end[,'r_end'],
         as.character = T
      )

      #--------- Repeat contexts ---------#
      if(verbose){ message("Calculating the number of copies of the indel sequence are present in the 3' flanking sequence...") }
      #' Calculate the number of copies of the indel sequence are present in the flank sequence
      #'
      #' @description Helper function for extractSigsIndel(). Scans along the flanking sequence for
      #' copies of the indel sequence using a sliding window with a length equal to the length of
      #' the indel sequence. The sliding window jumps in increments equal to the indel length.
      #' Scanning stops prematurely when the indel sequence does not match the one in the sliding
      #' window. The indel sequence itself counts as one copy.
      #'
      #' Scanning needs only to be done from the 5' -> 3' direction as the reported POS of an indel
      #' in a repeat region is always the first position of the repeat region.
      #'
      #' @param indel.seq The indel sequence as a character
      #' @param flank.seq The flanking sequence as a character
      #'
      #' @return An integer stating the number of copies found in the flanking sequence
      nCopiesAlongFlank <- function(indel.seq, flank.seq){
         # indel.seq = "T"
         # flank.seq = "TTTTGCG"

         indel_length <- nchar(indel.seq)
         rail_seq <- paste0(indel.seq, flank.seq)

         count <- 0
         seq_window <- substr(rail_seq, 1, indel_length)
         while(indel.seq == seq_window){
            count <- count + 1
            seq_window <- substr(
               rail_seq,
               count*indel_length + 1,
               count*indel_length + indel_length
            )
         }

         return(count)
      }
      n_copies_along_flank <- unlist(Map(nCopiesAlongFlank, variants$indel_seq, r_flank, USE.NAMES=F))

      #--------- Microhomology contexts ---------#
      if(verbose){ message("Calculating the (max) number of bases that are homologous to the 5'/3' flanking sequence...") }
      #' Calculate the number of bases that are homologous to the 3' flanking sequence
      #'
      #' @description Helper function for extractSigsIndel(). Scans a maximum of 1 indel length from
      #' 5' to 3' in the flanking sequence to detect identical bases. Stops prematurely if a
      #' non-identical base is found.
      #'
      #' DSBs can be repaired using 3' microhomology, which can be present relative to either the +
      #' or - strand. Therefore, both left and right flanking sequences of the indel need to be
      #' considered. When calculating left flanking homology (i.e. homology in the 3' direction for
      #' the - strand), the reverse complement of the indel sequence and flanking sequence need to
      #' be taken. However, the reverse can be taken for the calculation to save computation.
      #'
      #' @param indel.seq The indel sequence as a character
      #' @param flank.seq The flanking sequence as a character
      #'
      #' @return An integer stating the number of bases in microhomology
      nBasesMH <- function(indel.seq, flank.seq, indel.len){
         #indel.sequence = "CTA"
         #flank.sequence = "C"

         indel_len <- nchar(indel.seq)
         indel.seq <- unlist(strsplit(indel.seq, ''))
         flank.seq <- unlist(strsplit(flank.seq, ''))[1:indel_len]

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

      n_bases_mh <- unlist(Map(function(indel_seq, l_flank, r_flank){
         mh_l <- nBasesMH(reverse(indel_seq), reverse(l_flank))
         mh_r <- nBasesMH(indel_seq, r_flank)

         max(mh_l,mh_r)

      }, variants$indel_seq, l_flank, r_flank, USE.NAMES=F))

      #--------- Assign repeat, microhomology, or no context ---------#
      if(verbose){ message('Determining indel contexts...') }
      context <- unlist(Map(function(n_copies_along_flank, n_bases_mh, indel_len){
         if (n_copies_along_flank >= 2){
            if(indel_len < 50){ context <-'rep' }
            else { context <- 'mh' }

         } else if(n_copies_along_flank >= 1 && n_bases_mh >= 2) {
            context <- 'mh'
         } else if(n_copies_along_flank >= 1 && n_bases_mh >= 1 && indel_len > 3 ) {
            context <- 'mh'

         } else {
            context <- 'none'
         }

         return(context)

      }, n_copies_along_flank, n_bases_mh, variants$indel_len))

      #--------- Gather components for counting final contexts/signatures ---------#
      if(verbose){ message('Counting indel context occurrences...') }
      ## Slightly redundant (could have just assigned components to a dataframe). But easier to debug
      sig_parts <- data.frame(
         indel_type = variants$indel_type,
         context,
         indel_len = variants$indel_len,
         n_copies_along_flank,
         n_bases_mh
      )

      ## Bin values larger than cap into one bin for indel_len and n_bases_mh
      sig_parts <- within(sig_parts,{
         indel_len[indel_len >= indel.len.cap] <- indel.len.cap
         n_bases_mh[n_bases_mh >= n.bases.mh.cap] <- n.bases.mh.cap
      })

      ## Count occurrences of each signature
      sig_occurrences <- table(with(sig_parts,{
         unlist(Map(function(indel_type,context,indel_len,n_copies_along_flank,n_bases_mh){
            if(context == 'mh'){
               paste(indel_type, context, 'bimh', n_bases_mh, sep = '.')
            } else {
               paste(indel_type, context, 'len', indel_len, sep = '.')
            }
         },
         indel_type, context, indel_len, n_copies_along_flank, n_bases_mh))
      }))

      ## Fill in indel signature matrix that was initiated at the start of the function
      indel_sigs[names(sig_occurrences)] <- sig_occurrences
   }

   if(verbose){ message('Returning indel context counts...') }
   out <- matrix(indel_sigs, ncol = 1)
   rownames(out) <- names(indel_sigs)
   colnames(out) <- if(is.null(sample.name)){ basename(vcf.file) } else { sample.name }

   return(out)
}


