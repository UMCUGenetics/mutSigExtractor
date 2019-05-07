#' Extract relevant variant info for extracting SNV, indel, and SV signatures.
#'
#' @param vcf.file Path to the vcf file
#' @param mode A character stating which type of signature is to be extracted: 'snv','indel', or 'sv'
#' @param get.other.indel.allele Only applies when mode=='indel' For indels, some vcfs only report
#' the sequence of one allele (REF for deletions and ALT for insertions). If TRUE, the unreported
#' allele will be retrieved from the genome: a 5' base relative to the indel sequence. This base
#' will also be added to the indel sequence and the POS will be adjusted accordingly (POS=POS-1).
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
   vcf.file, mode=NULL, get.other.indel.allele=F, sv.caller='manta',
   ref.genome = DEFAULT_GENOME, chrom.group='all', vcf.filter=NULL, verbose=F
){

   if(!(mode %in% c('snv','indel','sv'))){ stop("Mode must be 'snv','indel', or 'sv'") }

   if(verbose){ message('Reading in vcf file...') }
   #vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/ICGC/vcf/OV-AU/snv_indel/AOCS-057_snv_indel.vcf.gz'
   #vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/ICGC/vcf/OV-AU/snv_indel/AOCS-001_snv_indel.vcf.gz'
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
   if(mode=='snv' || mode=='indel'){
      vcf_rr <- rowRanges(vcf)

      ## Unselect rows with multiple ALT sequences
      if(!identical(which(lengths(vcf_rr$ALT) > 1), integer(0))){
         if(verbose){ message('Some rows have multiple ALT sequences. These will be removed.') }
         vcf_rr <- vcf_rr[which(lengths(vcf_rr$ALT) == 1)]
      }

      ## Get main columns
      df <- data.frame(
         chrom=as.character(seqnames(vcf_rr)),
         pos=start(vcf_rr),
         ref=as.character(vcf_rr$REF),
         alt=as.character(unlist(vcf_rr$ALT)),
         stringsAsFactors=F
      )
   }

   #--------- SNV ---------#
   if(mode=='snv'){
      if(verbose){ message('Returning SNV trinucleotide contexts...') }
      snv_contexts <- data.frame(
         substitution = paste0(df$ref,'>',df$alt),
         tri_context = getSeq(
            x = eval(parse(text=ref.genome)),
            names = seqnames(vcf_rr),
            start = start(vcf_rr) - 1,
            end = end(vcf_rr) + 1,
            as.character = T
         ),
         stringsAsFactors = F
      )
      snv_contexts <- snv_contexts[nchar(df$ref)==1 & nchar(df$alt)==1,]
      return(snv_contexts)
   }

   #--------- Indel ---------#
   if(mode=='indel'){

      if(verbose){ message('Determining indel type...') }
      ## Calc sequence lengths
      df$ref_len <- nchar(df$ref)
      df$alt_len <- nchar(df$alt)

      ## Remove snvs
      df <- df[df$ref_len!=1 & df$alt_len!=1,]

      ## Determine indel type
      df$indel_type <- with(df,{
         unlist(Map(function(ref_len, alt_len){
            if(ref_len >= 2 & alt_len >= 2){
               if(ref_len == alt_len){ 'mnv_neutral' }
               else if(ref_len > alt_len){ 'mnv_del' }
               else if(ref_len < alt_len){ 'mnv_ins' }
            } else if(ref_len > alt_len){
               'del'
            } else if(ref_len < alt_len){
               'ins'
            }
         },ref_len, alt_len, USE.NAMES=F))
      })

      if(get.other.indel.allele){
         if(verbose){ message('Retrieving other indel allele...') }
         df_split <- lapply(
            list(del_type=c('del','mnv_del'),ins_type=c('ins','mnv_ins'),mnv_neutral='mnv_neutral'),
            function(i){ df[df$indel_type %in% i, ] }
         )

         ## Deletions
         ## ref:   'AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
         ## alt:  'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
         ## nchar('AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC') = nchar(ref) = 46
         ## getSeq(x=eval(parse(text = ref.genome)), names='chr4',start=84726292-1,84726292+46-1)
         ##       'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
         ## 5' base relative to ref ->
         ##    alt column
         ##    ref sequence
         df_split$del_type$alt <- with(df_split$del_type, {
            getSeq(
               x=eval(parse(text=ref.genome)),
               names=chrom, start=pos-1,end=pos-1,
               as.character=T
            )
         })
         df_split$del_type$ref <- with(df_split$del_type, { paste0(alt, ref) })
         df_split$del_type$pos <- df_split$del_type$pos-1

         ## Insertions
         ## ref:  ''
         ## alt:  'AGAGAGAGAGACAGAA'
         ## nchar('AGAGAGAGAGACAGAA') = nchar(alt) = 16
         ## getSeq(x=eval(parse(text = ref.genome)), names='chr12',start=6902128-1,6902128+16-1)
         ##      'GAGAGAGAGAGACAGAA'
         ## 5' base relative to alt ->
         ##    ref column
         ##    alt sequence
         ## Substract 1 from pos
         df_split$ins_type$ref <- with(df_split$ins_type, {
            getSeq(
               x=eval(parse(text=ref.genome)),
               names=chrom, start=pos-1,end=pos-1,
               as.character=T
            )
         })
         df_split$ins_type$alt <- with(df_split$ins_type, { paste0(ref,alt) })
         df_split$ins_type$pos <- df_split$ins_type$pos-1

         ## Unsplit df
         df <- do.call(rbind, df_split)
         rownames(df) <- NULL

         ## Recalculate ref/alt length
         df$ref_len <- nchar(df$ref)
         df$alt_len <- nchar(df$alt)
      }

      if(verbose){ message('Determining indel length and sequence...') }
      ## Determine indel length
      df$indel_len <- abs(df$alt_len-df$ref_len)

      ## Determine indel seq
      df$indel_seq <- with(df,{
         unlist(Map(function(ref,alt,indel_type,indel_len){
            indel_start_pos <- 2
            if(indel_type %in% c('del','mnv_del')){ ## dels
               substring(ref, indel_start_pos, indel_start_pos+indel_len-1)
            } else if(indel_type %in% c('ins','mnv_ins')){ ## ins
               substring(alt, indel_start_pos, indel_start_pos+indel_len-1)
            } else {
               NA
            }
         },ref,alt,indel_type,indel_len, USE.NAMES=F))
      })

      ## Output
      if(verbose){ message('Returning indel characteristics...') }
      out <- df[df$indel_type %in% c('ins','del'),]
      out <- out[,c('chrom','pos','ref','alt','indel_len','indel_type','indel_seq')]
      return(out)
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
         vcf_no_partners <- vcf[grepl('o', names(vcf)) | grepl('b', names(vcf))]

         df <- data.frame(
            id = names(vcf_no_partners),
            partner_type = grep(names(vcf_no_partners), '[ob]$', value=T),
            chrom_ref = gsub('chr','',as.character(seqnames(vcf_no_partners))),
            pos_ref = start(vcf_no_partners),
            seq_ref = as.character(rowRanges(vcf_no_partners)$REF),
            alt = as.character(rowRanges(vcf_no_partners)$ALT),
            stringsAsFactors = F
         )

         if(verbose){ message('Formatting ALT and calculating preliminary SV length...') }
         alt_split <- grep('[\\d\\w]+:\\d+',df$alt, value=T)
         alt_coord <- as.data.frame(do.call(lapply(alt_split, function(i){
            #i=alt_split[[1]]
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
