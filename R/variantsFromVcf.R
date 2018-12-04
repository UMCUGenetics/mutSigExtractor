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
   #vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_PCAGW/scripts/mutSigExtractor/R_source/gridss_output/XXXXXXXX.purple.sv.vcf.gz'
   vcf <- readVcf(vcf.file)

   ## Deal with empty vcfs
   if(identical(width(vcf),integer(0))){
      out <- NA

   } else {
      if(!is.null(ref.genome)){
         ref_genome <- eval(parse(text = ref.genome))
         ref_organism <- GenomeInfoDb::organism(ref_genome)
         seqlevelsStyle(vcf) <- seqlevelsStyle(ref_genome)
      }

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
         if (length(ref) == 0){ ref <- as.character(vcf$ref) } ## Allow both uppercase and lowercase column names when retrieving ref/alt

         alt <- as.character(unlist(vcf$ALT))
         if (length(alt) == 0){ alt <- vcf$alt }

         ## Unspecified mode
         if(is.null(mode)){
            if(verbose){ message('No mode specified. Outputting chrom, pos, ref, alt columns.') }
            out <- data.frame(
               chrom = chrom,
               pos = pos,
               ref = ref,
               alt = alt,
               #filter = vcf$FILTER,
               stringsAsFactors = F
            )
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

            ## Assign preliminary variant types and remove SNVs and MNVs
            ref_lengths <- nchar(ref)
            alt_lengths <- nchar(alt)
            indel_lengths <- abs(alt_lengths - ref_lengths)

            indel_type_sequence <- lapply(1:length(ref), function(i){
               ## Assign variant type
               if(ref_lengths[i] == 1 & alt_lengths[i] == 1){
                  variant_type <- 'snv'
               } else if(ref_lengths[i] >= 2 & alt_lengths[i] >= 2){
                  if(ref_lengths[i] == alt_lengths[i]){ variant_type <- 'mnv_neutral' }
                  else if(ref_lengths[i] > alt_lengths[i]){ variant_type <- 'mnv_del' }
                  else if(ref_lengths[i] < alt_lengths[i]){ variant_type <- 'mnv_ins' }
               } else if(ref_lengths[i] > alt_lengths[i]){
                  variant_type <- 'del'
               } else if(ref_lengths[i] < alt_lengths[i]){
                  variant_type <- 'ins'
               }

               ## Get indel seq
               indel_start_pos <- 2
               if(variant_type == 'del' ||  variant_type == 'mnv_del'){ ## dels/mnv
                  indel_seq <- substring(ref[i],indel_start_pos, indel_start_pos+indel_lengths[i]-1)
               } else if(variant_type == 'ins' ||  variant_type == 'mnv_ins'){ ## ins/mnv
                  indel_seq <- substring(alt[i],indel_start_pos, indel_start_pos+indel_lengths[i]-1)
               } else {
                  indel_seq <- NA
               }

               return( list(variant_type = variant_type, indel_seq = indel_seq) )
            })

            variant_type <- unlist(lapply(indel_type_sequence, function(i){ i$variant_type }))
            indel_seq <- unlist(lapply(indel_type_sequence, function(i){ i$indel_seq }))

            out <- data.frame(
               chrom = chrom,
               pos = pos,
               ref = ref,
               alt = alt,
               indel_len = indel_lengths,
               indel_type = variant_type,
               indel_seq = indel_seq,
               stringsAsFactors = F
            )

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

            vcf_as_df <- data.frame(
               id = names(vcf_no_partners),
               partner_type = unlist(str_extract_all(names(vcf_no_partners), '[ob]$')),
               chrom_ref = str_remove_all(as.character(seqnames(vcf_no_partners)), 'chr'),
               pos_ref = start(vcf_no_partners),
               seq_ref = as.character(rowRanges(vcf_no_partners)$REF),
               alt = as.character(rowRanges(vcf_no_partners)$ALT),
               stringsAsFactors = F
            )

            alt_split <- str_extract_all(vcf_as_df$alt, '[\\d\\w]+:\\d+')
            alt_coord <- lapply(alt_split, function(i){
               #i=alt_split[[1]]
               if(length(i) == 0){
                  c(NA,NA)
               } else {
                  unlist(str_split(i, ':'))
               }
            })
            alt_coord <- as.data.frame(do.call(rbind,alt_coord))
            colnames(alt_coord) <- c('chrom_alt','pos_alt')


            vcf_as_df <- cbind(vcf_as_df, alt_coord)
            vcf_as_df$pos_alt <- as.numeric(as.character(vcf_as_df$pos_alt))

            vcf_as_df$sv_len_pre <- vcf_as_df$pos_alt - vcf_as_df$pos_ref

            out <- do.call(rbind,lapply(1:nrow(vcf_as_df), function(i){
               row <- vcf_as_df[i,]

               if(row$partner_type == 'b'){
                  sv_type <- 'SGL'
               } else if(row$chrom_ref != row$chrom_alt){
                  sv_type <- 'TRA'
               } else if(row$sv_len_pre == 1){
                  sv_type <- 'INS'
               } else if(grepl('\\w+\\[.+\\[', row$alt)){
                  sv_type <- 'DEL'
               } else if(grepl('\\].+\\]\\w+', row$alt)){
                  sv_type <- 'DUP'
               } else if(grepl('\\w+\\].+\\]', row$alt) | grepl('\\[.+\\[\\w+', row$alt) ){
                  sv_type <- 'INV'
               } else {
                  sv_type <- NA
               }

               if(sv_type %in% c('SGL','TRA')){ sv_len <- NA }
               else{ sv_len <- row$sv_len_pre }

               return(data.frame(sv_type, sv_len, stringsAsFactors = F))
            }))

         } else {
            stop('Please specify SV caller')
         }

      }
   }

   return(out)
}
