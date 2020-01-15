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
   vcf.file, mode=NULL, sv.caller='gridss', ref.genome=DEFAULT_GENOME, chrom.group='all',
   vcf.filter=NA, verbose=F
){

   # mode='sv'
   # sv.caller='manta'
   # vcf.filter='PASS'
   # verbose=T
   # vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/example/vcf/PD3905_sv.vcf.gz'

   if(!(mode %in% c('snv_indel','sv'))){ stop("Mode must be 'snv_indel', or 'sv'") }

   if(verbose){ message('Reading in vcf file...') }

   vcf_fields <- c('CHROM','POS','REF','ALT','FILTER')
   if(mode=='sv'){
      if(sv.caller=='gridss'){
         vcf_fields <- c(vcf_fields, 'ID')
      } else if(sv.caller=='manta'){
         vcf_fields <- c(vcf_fields, 'INFO')
      }
   }

   vcf <- readVcfFields(vcf.file, vcf_fields)
   colnames(vcf) <- tolower(colnames(vcf))

   if(nrow(vcf)==0){
      if(verbose){ warning('VCF contains no rows. Returning NA') }
      return(NA)
   }

   ## Set chromosome names to the same used in the supplied ref genome
   vcf$chrom <- as.character(vcf$chrom)
   if(!is.null(ref.genome)){
      ref_genome <- eval(parse(text=ref.genome))
      ref_organism <- GenomeInfoDb::organism(ref_genome)
      seqlevelsStyle(vcf$chrom) <- seqlevelsStyle(ref_genome)
   }

   ## Keep certain chromosome types
   if(chrom.group != 'all'){
      if(is.null(ref_genome)){ stop('chromosome.group was specified but no reference genome was provided') }

      genome_chrom_group_names <- extractSeqlevelsByGroup(
         species=ref_organism,
         style=seqlevelsStyle(ref_genome),
         group=chrom.group
      )
      target_chrom_group_names <- intersect(genome_chrom_group_names, vcf$chrom)
      vcf$chrom <- vcf$chrom[vcf$chrom %in% target_chrom_group_names]
   }

   ## Filter vcf
   if(!is.na(vcf.filter)){
      if(verbose){ message('Only keeping variants where FILTER is ', paste(vcf.filter,collapse=', ')) }
      vcf <- vcf[vcf$filter %in% vcf.filter,]
   }

   if(nrow(vcf)==0){
      if(verbose){ warning('After filtering, VCF contains no rows. Returning NA') }
      return(NA)
   }

   #========= SNV/Indels =========#
   if(mode=='snv_indel'){

      ## Unselect rows with multiple ALT sequences
      if(verbose){ message('Removing rows with multiple ALT sequences...') }
      vcf <- vcf[!grepl(',',vcf$alt),]

      ## Post-processing
      vcf$filter <- NULL

      return(vcf)
   }

   #========= SV =========#
   if(!(sv.caller %in% c('manta','gridss'))){
      stop("Please specify valid SV caller: 'manta','gridss'")
   }

   if(sv.caller=='manta'){

      if(verbose){ message('Returning SV length and type...') }

      out <- getInfoValues(vcf$info,c('SVTYPE','SVLEN'))
      colnames(out) <- c('sv_type','sv_len')
      out$sv_len <- as.numeric(out$sv_len)
      out$sv_len[out$sv_type=='TRA'] <- NA

      return(out)
   }

   if(sv.caller=='gridss'){

      ## ID nomenclature:
      ## ends with o: 5' breakend
      ## ends with h: 3' breakend
      ## ends with b: unpaired breakend
      ## 'o' breakends have positive SV length

      if(verbose){ message('Keeping one breakend...') }
      df <- vcf[grepl('o$', vcf$id) | grepl('b$', vcf$id),]

      if(verbose){ message('Determining partner type...') }
      df$partner_type <- unlist(regmatches(df$id, gregexpr('[ob]$', df$id)))

      df$chrom_ref <- gsub('chr','',df$chrom)
      df$pos_ref <- df$pos

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
