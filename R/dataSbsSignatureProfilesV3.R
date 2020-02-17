#' SBS signature profiles (new; v3, May 2019)
#'
#' A matrix
#' rows: 96-trinucleotide context
#' cols: single base substitution (SBS) signatures
#'
#' Source: https://www.synapse.org/#!Synapse:syn12009743
#'
#' @docType data
#'
#' @usage data(SBS_SIGNATURE_PROFILES_V2)
'SBS_SIGNATURE_PROFILES_V2'

# ## Code to create RData
# sbs_ref <- read.csv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/raw/sigProfiler_SBS_signatures_2019_05_22.csv')
#
# rownames(sbs_ref) <- with(sbs_ref,{
#    paste0(
#       substr(SubType,1,1),
#       '[',Type,']',
#       substr(SubType,3,3)
#    )
# })
# sbs_ref <- sbs_ref[,grep('^SBS',colnames(sbs_ref))]
# sbs_ref <- as.matrix(sbs_ref)
#
# dim(sbs_ref)
#
# ## Remove signatures potentially related to sequencing artefacts
# excl_sbs_sigs <- c(
#    'SBS27','SBS43','SBS45','SBS46','SBS47',
#    'SBS48','SBS49','SBS50','SBS51','SBS52',
#    'SBS53','SBS54','SBS55','SBS56','SBS57',
#    'SBS58','SBS59','SBS60'
# )
#
# sbs_ref <- sbs_ref[,!(colnames(sbs_ref) %in% excl_sbs_sigs)]
#
# SBS_SIGNATURE_PROFILES_V3 <- sbs_ref
# save(SBS_SIGNATURE_PROFILES_V3, file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/SBS_SIGNATURE_PROFILES_V3.RData')
