use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.5'
.libPaths(use_library)

# BiocManager::install("EpiDISH")
library(EpiDISH)

# === INFO =====================================================================
# Cohort specific and mega matrix

data_dir <- '~/methylCHART/data'

dnam <- bigmemory::attach.big.matrix(
  file.path(data_dir, 'genr_dnam.desc'), readonly = TRUE, lockfile = TRUE)

which(!rownames(centUniLIFE.m) %in% colnames(dnam))

dnam <- dnam[, rownames(centUniLIFE.m)]


epid_funct <- epidish(beta.m = t(dnam[]), ref.m = centUniLIFE.m, method = "RPC")
# epid_unilife_funct <- epid_funct$estF              #get the estimated franctions
# saveRDS(epid_unilife_funct, "ell_counts_unilife_functional.rds")
