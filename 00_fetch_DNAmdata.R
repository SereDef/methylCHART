# library(haven)
# library(dplyr)
# library(jsonlite)

dnam_dir <- '~/GENR3/Methylation/Release4'

normalization <- "Functional"

genr <- haven::read_sav(file.path('data', 
                                  'DNAm_selection_file.sav'))

# Select cpgs 
set.seed(73)
cluster_rep_cpgs <- readRDS("~/MPSR/DNAm-clustering/metadata/CpG_metadata.rds") |>
  dplyr::slice_sample(n = 1, by = p2_cluster) |>
  dplyr::transmute(
    cpg = as.character(cpg),
    label = purrr::map(paste("Cluster", p2_cluster), ~ list(label = .x, pheno = NULL))
  ) |>
  (\(x) stats::setNames(x$label, x$cpg))()

target_cpgs <- c(list(
  cg00528572 = list(label='Maternal aging', pheno=c('mom_age', 'mom_age_bin')),  # Alex
  cg05575921 = list(label='Maternal smoking', pheno=c('smoke')),  # holy grail
  cg11945228 = list(label='GPF (not at birth, yes in childhood)', pheno=c('cbcl_tot')), 
  cg10644885 = list(label='Asthma (not at birth, yes in childhood)', pheno=c('wheez')), 
  cg22708087 = list(label='Asthma (positive at birth, negative in childhood)', pheno=c('wheez')), 
  cg25520701 = list(label='ADHD (CREB5; yes birth, not childhood)', pheno=c('cbcl_adhd')), 
  cg01271805 = list(label='ADHD (ERC2; yes birth, not childhood)', pheno=c('cbcl_adhd')), 
  cg09989037 = list(label='ADHD (ST3GAL3; yes birth, not childhood)', pheno=c('cbcl_adhd'))
  ), cluster_rep_cpgs)

# Some clean-up
genr$Period <- factor(genr$Period,
                      levels = c("Birth", "Age5", "Age9", "Age13", "Age17"))

genr$Array <- genr$Batch
genr$Batch <- interaction(genr$Period, genr$Array, drop = TRUE)

# Remove bridges 
genr_nobridges <- genr |>
  dplyr::filter(!(Bridge == "Bridge_Sample" & Array == 'EPICv2'))

# Check bridge removal 
table(genr[c("Period", "Array")])
table(genr_nobridges[c("Period", "Array")])

# Define all Period / Array combinations to read
subsets <- unique(genr_nobridges[c("Period", "Array")])
subsets

# Copy and set up column for cpg data
data <- genr_nobridges
# data$cpg <- NA_real_

# is.character(data$Sample_ID) # TRUE
rm(genr, genr_nobridges)

dnam_data <- data[, 'Sample_ID']
dnam_data[names(target_cpgs)] <- NA_real_

for (i in seq_len(nrow(subsets))) {
  period_name <- ifelse(subsets$Period[i] == "Birth", "Birth", 
                        paste0(subsets$Period[i],"y"))
  
  array_name  <- ifelse(subsets$Array[i] == "450k", "450K", 
                        subsets$Array[i])
  
  cli::cli_rule("Period: {.strong {period_name}} | Array: {.strong {array_name}}")
  
  rds_path <- list.files(
    path = file.path(dnam_dir, period_name, array_name, normalization), 
    pattern = "ALL", full.names = TRUE)
  
  if (length(rds_path) != 1L) {
    cli::cli_alert_warning("No or multiple matching RDS file found, skipping subset.")
    next
  }
  
  cli::cli_progress_step("Reading data", spinner = TRUE)
  
  # cpg_data <- as.data.frame(data.table::fread(rds_path))[target_cpg, ]
  cpg_data <- readRDS(rds_path)
  
  available_cpgs <- intersect(names(target_cpgs), rownames(cpg_data))
  if (!length(available_cpgs)) {
    cli::cli_alert_warning("None of the target CpGs are in this file.")
    cli::cli_progress_done()
    next
  } else if (length(available_cpgs) < length(target_cpgs)) {
    missing_cpgs <- setdiff(names(target_cpgs), available_cpgs)
    cli::cli_alert_warning(
      "{missing_cpgs} missing in this file.")
    cli::cli_alert_info('tag{?s}: {purrr::map_chr(target_cpgs[missing_cpgs], "label")}')
  }
  
  cpg_data <- cpg_data[available_cpgs, ]
  
  cli::cli_progress_done()
  
  # Match ID order and check there are no overlaps or missing matches 
  idx <- match(dnam_data$Sample_ID, colnames(cpg_data))
  
  has_match <- !is.na(idx)
  
  tot_samples <- ncol(cpg_data)
  tot_matched <- sum(has_match)
  cli::cli_inform("Total n: {tot_samples} | {tot_matched} matched in selection file.")
  
  already_set <- !is.na(dnam_data[has_match, names(target_cpgs)])
  
  if (any(already_set)) {
    cli::cli_alert_warning("{already_set} ID{?s} match multiple sets.")
  }
  
  dnam_data[has_match, available_cpgs] <- t(cpg_data[, idx[has_match]])

}

rm(subsets, cpg_data, already_set, has_match, idx, 
   array_name, period_name, rds_path, tot_matched, tot_samples, i)


summary(dnam_data)

all_na <- sapply(dnam_data, function(x) all(is.na(x)))
not_cpg <- names(dnam_data)[all_na] # CpGs not found

dnam_data[not_cpg] <- NULL
target_cpgs[not_cpg] <- NULL

write.csv(dnam_data, file.path('data', 'DNAm_data.csv'))
jsonlite::write_json(target_cpgs, "cpg_targets.json", auto_unbox = TRUE, pretty = TRUE)

# ==== Covariates / phenotypes =================================================

covs <- haven::read_sav(file.path('data',
                                  'Family_GeneralData_Pregnancy_20251111.sav')) |>
  dplyr::transmute(IDM,
                   sex = haven::as_factor(GENDERPREG), # # 1 = boy 2 = girl
                   mom_age = AGE_M_Birth, # years
                   mom_age_bin = as.factor(ifelse(mom_age > 35, '> 35', '=< 35'))
                   ) 

msmk <- haven::read_sav(file.path('data',
                                  'Mother_Smoking_Pregnancy_20260216.sav')) |>
  dplyr::transmute(IDM, smoke = haven::as_factor(SMOKE_ALL))
                        # 1 = never smoked during pregnancy
                        # 2 = smoked until pregnancy was known
                        # 3 = continued smoking in pregnancy

adhd <- haven::read_sav(file.path('data',
                                  'Child_CBCL_Quest09Y_20250805.sav')) |>
  dplyr::transmute(IDC, 
                   cbcl_att = sum_att_9m, 
                   cbcl_adhd = sumD_adh_9m,
                   cbcl_int = sum_int_9m, 
                   cbcl_ext = sum_ext_9m,
                   cbcl_tot = cbcl_sum_9m)

asma <- haven::read_sav(file.path('data',
                                  'Child_WheezingPatterns_0-4Y_20251002.sav')) |>
  dplyr::transmute(IDC, wheez = wh4y)

cov_data <- Reduce(function(x, y) merge(x, y, by = "IDM", all.x = TRUE), 
                   list(data, covs, msmk))
cov_data <- Reduce(function(x, y) merge(x, y, by = "IDC", all.x = TRUE), 
               list(cov_data, adhd, asma))

rm(covs, msmk, adhd, asma)

write.csv(cov_data, file.path('data', 'cov_data.csv'))
