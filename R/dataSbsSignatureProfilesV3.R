#' SBS signature profiles (new; v3, May 2019)
#'
#' A matrix
#' rows: 96-trinucleotide context
#' cols: single base substitution (SBS) signatures
#'
#' Source: https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx
#'
#' @docType data
#'
#' @usage data(SBS_SIGNATURE_PROFILES_V3)
'SBS_SIGNATURE_PROFILES_V3'


# ## Code to create RData --------
# df <- openxlsx::read.xlsx(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/COSMIC_Mutational_Signatures_v3.1.xlsx',
#    sheet='SBS_GRCh37'
# )
#
# df$Type2 <- paste0(
#    substring(df$Subtype,1,1),
#    '[',
#    df$Type,
#    ']',
#    substring(df$Subtype,3,3)
# )
#
# df <- df[order(df$Type),]
#
# # identical(
# #    df$Type2,
# #    rownames(SBS_SIGNATURE_PROFILES_V3)
# # )
#
# m <- df[,grep('^SBS',colnames(df))]
# rownames(m) <- df$Type2
#
# excl_sigs <- read.delim(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/sigs_v3.1_exclusion.txt'
# )
# excl_sigs <- subset(excl_sigs, nchar(exclude_reason)!=0, sig_name, drop=T)
# m <- m[,!(colnames(m) %in% excl_sigs)]
#
# SBS_SIGNATURE_PROFILES_V3 <- m
# save(
#    SBS_SIGNATURE_PROFILES_V3,
#    file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/SBS_SIGNATURE_PROFILES_V3.RData'
# )
