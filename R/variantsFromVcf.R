#' Extract relevant variant info for extracting SNV, indel, and SV signatures.
#'
#' @param vcf.file Path to the vcf file
#' @param mode A character stating which type of signature is to be extracted. SNV: 'snv', indel: 'indel', SV: 'sv'.
#' If unspecified, CHROM, POS, REF, and ALT will be returned
#' @param ref.genome A character naming the BSgenome reference genome. Default is "BSgenome.Hsapiens.UCSC.hg19". If another
#' reference genome is indicated, it will also need to be installed.
#' @param chrom.group chrom.group can be 'auto' for autosomes, 'sex' for sex chromosomes/allosomes, 'circular' for circular
#' chromosomes. The default is 'all' which returns all the chromosomes.
#' @param vcf.filter A character or character vector to specifying which variants to keep, corresponding to the values in
#' the vcf FILTER column
#' @param verbose Print progress messages?
#'
#' @return A data frame containing the relevant variant info for extracting the indicated signature type
#' @export
variantsFromVcf <- function(vcf.file, mode = NULL, sv.caller = 'manta', ref.genome = DEFAULT_GENOME,
                            chrom.group = 'all', vcf.filter = NULL, verbose = F){

   vcf <- readVcf(vcf.file)

   ## Deal with empty vcfs
   if( identical(width(vcf),integer(0)) ){
      warning('VCF contains no rows. Returning NA')
      stop(return(NA))
   }

   if(!is.null(ref.genome)){
      ref_genome <- eval(parse(text = ref.genome))
      ref_organism <- GenomeInfoDb::organism(ref_genome)
      seqlevelsStyle(vcf) <- seqlevelsStyle(ref_genome)
   }

   ## Remove non autosomal chromosomes?
   if(chrom.group != 'all'){
      if(is.null(ref_genome)){ stop('chromosome.group was specified but no reference genome was provided') }

      genome_chrom_group_names <- extractSeqlevelsByGroup(
         species = ref_organism,
         style = seqlevelsStyle(ref_genome),
         group = chrom.group)
      target_chrom_group_names <- intersect(genome_chrom_group_names, seqlevels(vcf))

      vcf <- keepSeqlevels(vcf, target_chrom_group_names, pruning.mode = "coarse")
   }

   ## Filter vcf
   if(!is.null(vcf.filter)){
      if(verbose){ message('Only keeping variants where FILTER %in% c(', paste(vcf.filter,collapse=', '), ')' ) }
      vcf <- vcf[which(rowRanges(vcf)$FILTER %in% vcf.filter)]
   }

   if( identical(width(vcf),integer(0)) ){
      warning('After filtering, VCF contains no rows. Returning NA')
      stop(return(NA))
   }

   #========= SNV/Indels =========#
   if(mode == 'snv' || mode == 'indel'){
      vcf <- rowRanges(vcf)

      ## Unselect rows with multiple ALT sequences
      if(!identical(which(lengths(vcf$ALT) > 1), integer(0))){
         if(verbose){ message('Some rows have multiple ALT sequences. These will be removed.') }
         vcf <- vcf[which(lengths(vcf$ALT) == 1)]
      }

      ## Get main columns
      chrom <- as.character(seqnames(vcf))
      pos <- start(vcf)

      ref <- as.character(vcf$REF)
      alt <- as.character(unlist(vcf$ALT))

      ## Unspecified mode
      if(is.null(mode)){
         if(verbose){ message('No mode specified. Outputting chrom, pos, ref, alt columns.') }
         out <- data.frame(chrom,pos,ref,alt,stringsAsFactors = F)
      }

      #--------- SNV ---------#
      else if(mode == 'snv'){
         out <- data.frame(
            substitution = paste0(ref,'>',alt),
            tri_context = getSeq(
               x = eval(parse(text=ref.genome)),
               names = seqnames(vcf),
               start = start(vcf) - 1,
               end = end(vcf) + 1,
               as.character = T
            ),
            stringsAsFactors = F
         )
         out <- out[nchar(ref) == 1 & nchar(alt) == 1,]
      }

      #--------- Indel ---------#
      else if(mode == 'indel'){
         indel_type_sequence <-  Map(function(ref, alt, ref_lengths, alt_lengths, indel_lengths){
            ## Assign variant type
            if(ref_lengths == 1 & alt_lengths == 1){
               variant_type <- 'snv'
            } else if(ref_lengths >= 2 & alt_lengths >= 2){
               if(ref_lengths == alt_lengths){ variant_type <- 'mnv_neutral' }
               else if(ref_lengths > alt_lengths){ variant_type <- 'mnv_del' }
               else if(ref_lengths < alt_lengths){ variant_type <- 'mnv_ins' }
            } else if(ref_lengths > alt_lengths){
               variant_type <- 'del'
            } else if(ref_lengths < alt_lengths){
               variant_type <- 'ins'
            }

            ## Get indel seq
            indel_start_pos <- 2
            if(variant_type == 'del' ||  variant_type == 'mnv_del'){ ## dels/mnv
               indel_seq <- substring(ref,indel_start_pos, indel_start_pos+indel_lengths-1)
            } else if(variant_type == 'ins' ||  variant_type == 'mnv_ins'){ ## ins/mnv
               indel_seq <- substring(alt,indel_start_pos, indel_start_pos+indel_lengths-1)
            } else {
               indel_seq <- NA
            }

            ## Return a vector instead of a list to prevent dataframe having lists as columns
            ## indel_lengths will be a character, but is converted into an integer below
            return( data.frame(indel_lengths, variant_type, indel_seq) )
         },
            ref, alt, ref_lengths = nchar(ref), alt_lengths = nchar(alt),
            indel_lengths = abs(nchar(alt)-nchar(ref)), USE.NAMES = F
         )
         indel_type_sequence <- as.data.frame(do.call(rbind, indel_type_sequence), stringsAsFactors=F)
         colnames(indel_type_sequence) <- c('indel_len','indel_type','indel_seq')
         indel_type_sequence$indel_seq <- as.character(indel_type_sequence$indel_seq)

         out <- cbind(chrom, pos, ref, alt, indel_type_sequence, stringsAsFactors = F)
         #out <- out[!(out$indel_type %in% c('snv','mnv_neutral')),]
         out <- out[out$indel_type %in% c('del','ins'),]
      }
   }

   #========= SV =========#
   else if (mode == 'sv'){

      if(sv.caller == 'manta'){
         sv_len <- unlist(lapply(info(vcf)$SVLEN, function(i){
            if(length(i) == 0){ i <- NA }
            else { i }
         }))
         sv_len <- abs(sv_len) ## Convert negative sv_len from DEL to positive

         out <- data.frame(
            sv_type = info(vcf)$SVTYPE,
            sv_len = sv_len
         )

      } else if(sv.caller == 'gridss'){
         ## Retrieve sense partners and unpartnered variants
         vcf_no_partners <- vcf[grepl('o', names(vcf)) | grepl('b', names(vcf))]

         df <- data.frame(
            id = names(vcf_no_partners),
            partner_type = unlist(str_extract_all(names(vcf_no_partners), '[ob]$')),
            chrom_ref = str_remove_all(as.character(seqnames(vcf_no_partners)), 'chr'),
            pos_ref = start(vcf_no_partners),
            seq_ref = as.character(rowRanges(vcf_no_partners)$REF),
            alt = as.character(rowRanges(vcf_no_partners)$ALT),
            stringsAsFactors = F
         )

         ## Preprocessing before decision tree
         alt_split <- str_extract_all(df$alt, '[\\d\\w]+:\\d+')
         alt_coord <- lapply(alt_split, function(i){
            #i=alt_split[[1]]
            if(length(i) == 0){ c(NA,NA) }
            else { unlist(str_split(i, ':')) }
         })
         alt_coord <- as.data.frame(do.call(rbind,alt_coord))
         colnames(alt_coord) <- c('chrom_alt','pos_alt')

         df <- cbind(df, alt_coord)
         df$pos_alt <- as.numeric(as.character(df$pos_alt))

         df$sv_len_pre <- df$pos_alt - df$pos_ref

         ## Decision tree
         out <- do.call(rbind,Map(function(partner_type, chrom_ref, chrom_alt, alt, sv_len_pre){
            if(partner_type == 'b'){
               sv_type <- 'SGL'
            } else if(chrom_ref != chrom_alt){
               sv_type <- 'TRA'
            } else if(sv_len_pre == 1){
               sv_type <- 'INS'
            } else if(grepl('\\w+\\[.+\\[', alt)){
               sv_type <- 'DEL'
            } else if(grepl('\\].+\\]\\w+', alt)){
               sv_type <- 'DUP'
            } else if(grepl('\\w+\\].+\\]', alt) | grepl('\\[.+\\[\\w+', alt) ){
               sv_type <- 'INV'
            } else {
               sv_type <- NA
            }

            if(sv_type %in% c('SGL','TRA')){ sv_len <- NA }
            else{ sv_len <- sv_len_pre }

            return(data.frame(sv_type, sv_len, stringsAsFactors = F))
         }, df$partner_type, df$chrom_ref, df$chrom_alt, df$alt, df$sv_len_pre, USE.NAMES = F))

      } else {
         stop('Please specify SV caller')
      }

   }

   return(out)
}
