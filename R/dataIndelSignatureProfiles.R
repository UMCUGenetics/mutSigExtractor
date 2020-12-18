#' Indel signature profiles
#'
#' A matrix
#' rows: indel context
#' cols: indel signature
#'
#' Source: https://cancer.sanger.ac.uk/cosmic/signatures/ID/index.tt
#' Each signature profile was downloaded manually, e.g. at https://cancer.sanger.ac.uk/sigs-assets-20/ID_vignettes/sigProfiler_ID_signatures_ID1.csv
#'
#' @docType data
#'
#' @usage data(INDEL_SIGNATURE_PROFILES)
'INDEL_SIGNATURE_PROFILES'

# ## Code to create RData --------
# ## Gather CSV files
# csv_files <- list.files(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/raw/indel_signatures',
#    pattern='sigProfiler.*.csv$',full.names=T
# )
# csv_files <- naturalsort::naturalsort(csv_files)
#
# l_sig_profiles <- lapply(csv_files, read.csv, stringsAsFactors=F)
#
# m <- do.call(cbind, lapply(l_sig_profiles,`[`,2))
# colnames(m) <- gsub('_\\w+$','',colnames(m))
# rownames(m) <- l_sig_profiles[[1]]$MutationType
# m <- m[grep(':',rownames(m)),]
#
# write.table(
#    m,
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/raw/sigProfiler_ID_signatures_3.1.txt',
#    sep='\t', quote=F
# )
#
# ## !!! Added mut_type_2 by hand
# m2 <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/raw/sigProfiler_ID_signatures_3.1.txt')
# m2 <- m2[!is.na(m2$mut_type_2),]
# context_names <- m2$mut_type_2
#
# m2 <- m2[,grep('^ID',colnames(m2))]
# rownames(m2) <- context_names
#
# INDEL_SIGNATURE_PROFILES <- as.matrix(m2)
# save(
#    INDEL_SIGNATURE_PROFILES,
#    file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/INDEL_SIGNATURE_PROFILES.RData'
# )
