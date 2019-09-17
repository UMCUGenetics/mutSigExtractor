#' Extract relevant variant info for extracting SNV, indel, and SV signatures.
#'
#' @param vcf.file Path to the vcf file
#' @param mode A character stating which type of signature is to be extracted: 'snv','indel', or 'sv'
#' @param sv.caller Only applies when mode=='sv'. At this moment supports 'manta' or 'gridss'.
#' Currently there is no standard to how SVs are reported in vcfs. Therefore, output from different
#' callers will need to be parsed separately.
#' @param ref.genome A character naming the BSgenome reference genome. Default is
#' "BSgenome.Hsapiens.UCSC.hg19". If another reference genome is indicated, it will also need to be
#' installed.
#' @param chrom.group chrom.group can be 'auto' for autosomes, 'sex' for sex chromosomes/allosomes,
#' 'circular' for circular chromosomes. The default is 'all' which returns all the chromosomes.
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param verbose Print progress messages?
#'
#' @return A data frame containing the relevant variant info for extracting the indicated signature type
#' @export
variantsFromVcf <- function(
   vcf.file, mode=NULL, sv.caller='manta', ref.genome=DEFAULT_GENOME, chrom.group='all',
   vcf.filter=NULL, verbose=F
){

   if(!(mode %in% c('snv_indel','sv'))){ stop("Mode must be 'snv_indel', or 'sv'") }

   if(verbose){ message('Reading in vcf file...') }

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
   if(mode=='snv_indel'){
      vcf_rr <- rowRanges(vcf)

      ## Unselect rows with multiple ALT sequences
      if(!identical(which(lengths(vcf_rr$ALT) > 1), integer(0))){
         if(verbose){ message('Removing rows with multiple ALT sequences...') }
         vcf_rr <- vcf_rr[which(lengths(vcf_rr$ALT) == 1)]
      }

      ## Get main columns
      bed <- data.frame(
         chrom=as.character(seqnames(vcf_rr)),
         pos=start(vcf_rr),
         ref=as.character(vcf_rr$REF),
         alt=as.character(unlist(vcf_rr$ALT)),
         stringsAsFactors=F
      )

      return(bed)
   }

   #========= SV =========#
   if (mode=='sv'){

      if(!(sv.caller %in% c('manta','gridss'))){
         stop("Please specify valid SV caller: 'manta','gridss'")
      }

      if(sv.caller=='manta'){
         if(verbose){ message('Returning SV length and type...') }
         sv_len <- unlist(lapply(info(vcf)$SVLEN, function(i){
            if(length(i) == 0){ i <- NA }
            else { i }
         }))
         sv_len <- abs(sv_len) ## Convert negative sv_len from DEL to positive

         out <- data.frame(
            sv_type = info(vcf)$SVTYPE,
            sv_len = sv_len
         )
         return(out)
      }

      if(sv.caller=='gridss'){
         if(verbose){ message('Retrieving sense partners and unpartnered SVs...') }
         vcf_no_partners <- vcf[grepl('o$', names(vcf)) | grepl('b$', names(vcf))]

         df <- data.frame(
            id = names(vcf_no_partners),
            partner_type = unlist(regmatches(names(vcf_no_partners), gregexpr('[ob]$', names(vcf_no_partners)))),
            chrom_ref = gsub('chr','',as.character(seqnames(vcf_no_partners))),
            pos_ref = start(vcf_no_partners),
            seq_ref = as.character(rowRanges(vcf_no_partners)$REF),
            alt = as.character(rowRanges(vcf_no_partners)$ALT),
            stringsAsFactors = F
         )

         if(verbose){ message('Formatting ALT and calculating preliminary SV length...') }
         alt_coord <- regmatches(df$alt, gregexpr('\\d|\\w+:\\d+', df$alt))
         alt_coord <- as.data.frame(do.call(rbind, lapply(alt_coord, function(i){
            if(length(i) == 0){ c(NA,NA) }
            else { unlist(strsplit(i, ':')) }
         })))
         colnames(alt_coord) <- c('chrom_alt','pos_alt')

         df <- cbind(df, alt_coord)
         df$pos_alt <- as.numeric(as.character(df$pos_alt))

         df$sv_len_pre <- df$pos_alt - df$pos_ref

         ## Decision tree
         if(verbose){ message('Determining and returning SV length and type...') }
         out <- with(df,{
            do.call(rbind,Map(function(partner_type, chrom_ref, chrom_alt, alt, sv_len_pre){
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
            }, partner_type, chrom_ref, chrom_alt, alt, sv_len_pre, USE.NAMES = F))
         })

         return(out)
      }
   }
}
