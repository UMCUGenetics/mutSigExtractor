#' Extract relevant variant info for extracting SNV, indel, and SV signatures.
#'
#' @param vcf.file Path to the vcf file
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param keep.chroms A character vector specifying which chromosomes to keep (chromosome names
#' should be in the style of the vcf). To keep autosomal and sex chromosomes for example use:
#' keep.chroms=c(1:22,'X','Y')
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param vcf.fields A character vector specifying the vcf columns to retrieve
#' @param verbose Print progress messages?
#'
#' @return A data frame containing the relevant variant info for extracting the indicated signature type
#' @export
variantsFromVcf <- function(
   vcf.file, ref.genome=DEFAULT_GENOME, keep.chroms=NULL,
   vcf.filter=NA, vcf.fields=c('CHROM','POS','REF','ALT','FILTER'), verbose=F
){
   # mode='sv'
   # sv.caller='manta'
   # vcf.filter='PASS'
   # verbose=T
   # vcf.file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/example/vcf/PD3905_snv_indel.vcf.gz'

   if(verbose){ message('Reading in vcf file...') }
   vcf <- readVcfFields(vcf.file, vcf.fields)
   colnames(vcf) <- tolower(colnames(vcf))

   if(nrow(vcf)==0){
      if(verbose){ warning('VCF contains no rows. Returning NA') }
      return(NA)
   }

   ## Keep certain chromosome types
   if(!is.null(keep.chroms)){
      if(verbose){ 'Only keeping chromosomes as indicated in keep.chroms...' }
      ## Force chromosome name style to that in ref genome
      vcf <- vcf[vcf$chrom %in% keep.chroms,]
   }

   ## Set chromosome names to the same used in the supplied ref genome
   vcf$chrom <- as.character(vcf$chrom)
   if(!is.null(ref.genome)){
      if(verbose){ 'Converting chrom name style to style in ref.genome...' }
      GenomeInfoDb::seqlevelsStyle(vcf$chrom)<- GenomeInfoDb::seqlevelsStyle(ref.genome)
   }

   ## Filter vcf
   if(!is.na(vcf.filter)){
      if(verbose){ message('Only keeping variants where FILTER is ', paste(vcf.filter,collapse=', ')) }
      vcf <- vcf[vcf$filter %in% vcf.filter,]
      vcf$filter <- NULL
   }

   if(nrow(vcf)==0){
      if(verbose){ warning('After filtering, VCF contains no rows. Returning NA') }
      return(NA)
   }

   return(vcf)
}
