use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.5'
.libPaths(use_library)

library(bigmemory)
# library(foreach)
# library(dplyr)

# === INFO =====================================================================
# Save both cohort specific and mega matrix

data_dir <- '~/methylCHART/data'

# TMP only select a subset of CpGs, eventually smart chucking goes here --------

pheno <- read.csv(file.path(data_dir, 'pheno_data.csv'), row.names = 1)
cpg_list <- readLines(file.path(data_dir, 'cpg_list.txt'))

# == CPG selection =============================================================
# Birth cluster representaiton 
set.seed(73)
cluster_rep_cpgs <- readRDS("~/MPSR/DNAm-clustering/metadata/CpG_metadata.rds") |>
  dplyr::filter(cpg %in% cpg_list) |> 
  dplyr::slice_sample(n = 1, by = p2_cluster) |>
  dplyr::transmute(
    cpg = as.character(cpg),
    label = purrr::map(paste("Cluster", p2_cluster), ~ list(label = .x, pheno = NULL))
  ) |>
  (\(x) stats::setNames(x$label, x$cpg))()

target_cpgs <- c(list(
  cg00528572 = list(label='Maternal aging', pheno=c('mom_age', 'mom_age_bin')),  # Alex
  cg05575921 = list(label='Maternal smoking', pheno=c('smoke')),  # holy grail
  
  # cg11945228 = list(label='GPF (not at birth, yes in childhood)', pheno=c('cbcl_tot')), 
  # -- removing becasue not in the shared cpg list
  
  cg10644885 = list(label='Asthma (not at birth, yes in childhood)', pheno=c('wheez')), 
  cg22708087 = list(label='Asthma (positive at birth, negative in childhood)', pheno=c('wheez')), 
  cg25520701 = list(label='ADHD (CREB5; yes birth, not childhood)', pheno=c('cbcl_adhd')), 
  cg01271805 = list(label='ADHD (ERC2; yes birth, not childhood)', pheno=c('cbcl_adhd')), 
  cg09989037 = list(label='ADHD (ST3GAL3; yes birth, not childhood)', pheno=c('cbcl_adhd')),
  cg06536614 = list(label='ASDf_nc886 1', pheno=c('cbcl_tot')),
  cg04481923 = list(label='ASDf_nc886 2', pheno=c('cbcl_tot')),
  cg18678645 = list(label='ASDf_nc886 3', pheno=c('cbcl_tot')),
  cg25652859 = list(label='SCZc_GNAS', pheno=c('cbcl_tot')),
  cg03836171 = list(label='SCZm_C7orf50', pheno=c('cbcl_tot')),
  # cg21870668 = list(label='SCZf_NoMHC_SNRPN', pheno=c('cbcl_tot')), 
  # --- removing this becasue 3691 (12% of the sample) is missing -- resolve the problem later 
  cg23248424 = list(label='ADHDc_GFPT2', pheno=c('cbcl_adhd')) 
  ), cluster_rep_cpgs)

target_cpgs_present <- target_cpgs[names(target_cpgs) %in% cpg_list]

# == Read data subset =============================================================
genr_dnam <- bigmemory::attach.big.matrix(
  file.path(data_dir, 'genr_dnam.desc'), readonly = TRUE, lockfile = TRUE)

alsp_dnam <- bigmemory::attach.big.matrix(
  file.path(data_dir, 'alsp_dnam.desc'), readonly = TRUE, lockfile = TRUE)

genr_sub <- genr_dnam[, names(target_cpgs_present)]
alsp_sub <- alsp_dnam[, names(target_cpgs_present)]

dnam_data <- rbind(genr_sub, alsp_sub)

# what to do with missing data mmm
colSums(is.na(dnam_data))

# Is this a few subjects
table(rowSums(is.na(dnam_data)))

# yes... Mean impute for now, will handle this later
na_idx <- which(is.na(dnam_data), arr.ind = TRUE)
dnam_data[na_idx] <- colMeans(dnam_data, na.rm = TRUE)[na_idx[, 2]]

colSums(is.na(dnam_data))

# reorder data to match pheno --------------------------------------------------
idx <- match(pheno$Sample_ID, rownames(dnam_data))
stopifnot(!anyNA(idx))
dnam_data <- dnam_data[idx, ]

identical(rownames(dnam_data), pheno$Sample_ID)   # TRUE: rows are aligned

summary(dnam_data)

# Turn to dataframe (TMP) for correct ID matching
dnam_df <- data.frame(Sample_ID = rownames(dnam_data), dnam_data, check.names = FALSE, row.names = NULL)

identical(dnam_df$Sample_ID, pheno$Sample_ID)

write.csv(dnam_df, file.path(data_dir, 'dnam_data.csv'))

jsonlite::write_json(target_cpgs_present, "cpg_targets.json", 
                     auto_unbox = TRUE, pretty = TRUE)
