#' DBS signature profiles
#'
#' A matrix
#' rows: dinucleotide base substitution (DBS) context
#' cols: DBS signatures
#'
#' Source: https://www.synapse.org/#!Synapse:syn12009743
#'
#' @docType data
#'
#' @usage data(DBS_SIGNATURE_PROFILES)
'DBS_SIGNATURE_PROFILES'

## Code to create RData
# dbs_ref <- read.csv('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/raw/sigProfiler_DBS_signatures.csv')
#
# rownames(dbs_ref) <- dbs_ref$Mutation.Type
# dbs_ref$Mutation.Type <- NULL
#
# DBS_SIGNATURE_PROFILES <- as.matrix(dbs_ref)
# save(DBS_SIGNATURE_PROFILES, file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/DBS_SIGNATURE_PROFILES.RData')
