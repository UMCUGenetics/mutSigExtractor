#' Extract indel signatures
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
#' @param ... Other arguments that can be passed to variantsFromVcf()
#'
#' @return A 1-column matrix containing the absolute indel signature contributions (i.e. the number of mutations contributing
#' to each mutational signature)
#' @export

extractSigsIndel <- function(vcf.file, sample.name = NULL, ref.genome = DEFAULT_GENOME,
                             indel.len.cap = 5, n.bases.mh.cap = 5, ...){
   ## Read vcf file
   variants <- variantsFromVcf(vcf.file, mode = 'indel', ref.genome, ...)
   #variants <- variantsFromVcf(vcf.file, mode = 'indel', ref.genome, vcf.filter = 'PASS')

   ## Initiate indel sig matrix
   # indel_sigs_grouped <- list(
   #    del.rep = structure(rep(0,indel.len.cap), names = paste0('del.rep.len.', 1:indel.len.cap)),
   #    ins.rep = structure(rep(0,indel.len.cap), names = paste0('ins.rep.len.', 1:indel.len.cap)),
   #    del.mh = structure(rep(0,n.bases.mh.cap), names = paste0('del.mh.bimh.', 1:n.bases.mh.cap)),
   #    ins.mh = structure(rep(0,n.bases.mh.cap), names = paste0('ins.mh.bimh.', 1:n.bases.mh.cap)),
   #    del.none = structure(rep(0,indel.len.cap), names = paste0('del.none.len.', 1:indel.len.cap)),
   #    ins.none = structure(rep(0,indel.len.cap), names = paste0('ins.none.len.', 1:indel.len.cap))
   # )

   indel_sigs <- c(
      paste0('del.rep.len.', 1:indel.len.cap),
      paste0('ins.rep.len.', 1:indel.len.cap),
      paste0('del.mh.bimh.', 1:n.bases.mh.cap),
      paste0('ins.mh.bimh.', 1:n.bases.mh.cap),
      paste0('del.none.len.', 1:indel.len.cap),
      paste0('ins.none.len.', 1:indel.len.cap)
   )

   indel_sigs <- structure(rep(0,length(indel_sigs)), names = indel_sigs)

   ## Pre-calculations for repeat and microhomology contexts
   if(is.data.frame(variants) && !is.na(variants)){ ## Deal with empty vcfs

      ## Getting flank start/end positions first, then retrieving the sequences using getSeq with a vectorized input
      ## improves speed significantly
      flanks_start_end <- do.call(rbind, lapply(1:nrow(variants), function(i){
         indelSeqFlanksStartEnd(variants[i,'chrom'],
                                variants[i,'pos'],
                                variants[i,'indel_len'],
                                variants[i,'indel_type'],
                                ## Cap n.indel.lengths.r to 3. The max n_copies_along_flank condition used below caps at >=2.
                                ## This improves speed significantly
                                n.indel.lengths.r = 3)
      }))

      l_flank <- getSeq(x = eval(parse(text = ref.genome)),
                        names = variants$chrom,
                        start = flanks_start_end[,'l_start'],
                        end = flanks_start_end[,'l_end'],
                        as.character = T)

      r_flank <- getSeq(x = eval(parse(text = ref.genome)),
                        names = variants$chrom,
                        start = flanks_start_end[,'r_start'],
                        end = flanks_start_end[,'r_end'],
                        as.character = T)


      ## Repeat contexts:
      ## Scanning needs only to be done from the 5' -> 3' direction. In a vcf, the reported POS of an indel in a
      ## repeat region is always the first position of the repeat region.
      n_copies_along_flank <- unlist(lapply(1:nrow(variants), function(i){
         nCopiesAlongFlank(variants[i,'indel_seq'], r_flank[i])
      }))

      ## Microhomology contexts:
      ## DSBs can be repaired using 3' microhomology, which can be present relative to either the + or - strand.
      ## Therefore, both left and right flanking sequences of the indel need to be considered.
      ## When calculating left flanking homology (i.e. homology in the 3' direction for the - strand), the reverse
      ## complement of the indel sequence and flanking sequence need to be taken. However, the reverse can be taken
      ## for the calculation to save computation.
      n_bases_mh <- do.call(rbind, lapply(1:nrow(variants), function(i){
         mh_l <- nBasesMH(reverse(variants[i,'indel_seq']), reverse(l_flank[i]))
         mh_r <- nBasesMH(variants[i,'indel_seq'], r_flank[i])

         #list(n_bases_mh_l = mh_l, n_bases_mh_r = mh_r)
         max(mh_l,mh_r)
      }))

      ## Assign repeat, microhomology, or no context
      context <- unlist(lapply(1:nrow(variants), function(i){
         n_copies_along_flank <- n_copies_along_flank[i]
         n_bases_mh <- n_bases_mh[i]
         indel_len <- variants[i,'indel_len']

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
      }))

      ## Gather components for counting final signatures
      sig_components <- data.frame(
         indel_type = variants$indel_type,
         context,
         indel_len = variants$indel_len,
         n_copies_along_flank,
         n_bases_mh
      )

      # ## for debugging
      # n_bases_mh <- as.data.frame(n_bases_mh)
      # sig_components <- cbind(variants,
      #                         l_flank = l_flank,
      #                         r_flank = r_flank,
      #                         n_bases_mh_l = unlist(n_bases_mh$n_bases_mh_l),
      #                         n_bases_mh_r = unlist(n_bases_mh$n_bases_mh_r),
      #                         n_copies_along_flank,
      #                         context
      #                         )
      #
      # write.table(sig_components, '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mutSigExtractor/R_source/call_indel_signatures_v2/XXXXXXXX.txt', sep = '\t')

      ## Bin values larger than cap into one bin for indel_len and n_bases_mh
      sig_components$indel_len[sig_components$indel_len >= indel.len.cap] <- indel.len.cap
      sig_components$n_bases_mh[sig_components$n_bases_mh >= n.bases.mh.cap] <- n.bases.mh.cap

      ## Count occurrences of each signature
      sig_occurrences <- apply(sig_components,1,function(i){
         indel_type <- i[1]
         context <- i[2]
         indel_len <- i[3]
         n_copies_along_flank <- i[4]
         n_bases_mh <- i[5]

         if(context == 'mh'){
            paste(indel_type, context, 'bimh', n_bases_mh, sep = '.')
         } else {
            paste(indel_type, context, 'len', indel_len, sep = '.')
         }
      })

      sig_occurrences <- table(sig_occurrences)

      ## Make final indel signature matrix based on the one initiated at the start of the function (with all zeroes)
      non_zero_sigs <- names(indel_sigs)[names(indel_sigs) %in% names(sig_occurrences)]
      indel_sigs[non_zero_sigs] <- sig_occurrences[non_zero_sigs]
   }

   out <- matrix(indel_sigs, ncol = 1)
   rownames(out) <- names(indel_sigs)

   if(is.null(sample.name)){ colnames(out) <- basename(vcf.file) }
   else { colnames(out) <- sample.name }

   return(out)
}


