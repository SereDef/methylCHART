# library(haven)
# library(dplyr)

dnam_dir <- '~/GENR3/Methylation/Release4'

target_cpg <- "cg05575921" # holy grail 
normalization <- "Functional"

genr <- haven::read_sav(file.path('data', 
                                  'DNAm_selection_file.sav'))

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
data$cpg <- NA_real_

# is.character(data$Sample_ID) # TRUE
rm(genr, genr_nobridges)

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
  cpg_data <- readRDS(rds_path)[target_cpg, ]
  cli::cli_progress_done()
  
  # Match ID order and check there are no overlaps or missing matches 
  idx <- match(data$Sample_ID, names(cpg_data))
  has_match <- !is.na(idx)
  
  tot_samples <- length(cpg_data)
  tot_matched <- sum(has_match)
  cli::cli_inform("Total n: {tot_samples} | {tot_matched} matched in selection file.")
  
  already_set <- !is.na(data$cpg[has_match])
  
  if (any(already_set)) {
    cli::cli_alert_warning("{already_set} ID{?s} match multiple sets.")
  }
  
  # Write cpg value
  data$cpg[has_match] <- cpg_data[idx[has_match]]

}

rm(subsets, cpg_data, already_set, has_match, idx, 
   array_name, period_name, rds_path, tot_matched, tot_samples, i)

# Add covariate information
covs <- haven::read_sav(file.path('data',
                                  'Family_GeneralData_Pregnancy_20251111.sav')) |>
  dplyr::transmute(IDM, sex = haven::as_factor(GENDERPREG)) # 1 = boy 2 = girl

outc <- haven::read_sav(file.path('data',
                                  'Mother_Smoking_Pregnancy_20260216.sav')) |>
  dplyr::transmute(IDM, smoke = haven::as_factor(SMOKE_ALL))
                        # 1 = never smoked during pregnancy
                        # 2 = smoked until pregnancy was known
                        # 3 = continued smoking in pregnancy

data <- Reduce(function(x, y) merge(x, y, by = "IDM", all.x = TRUE), 
                   list(data, covs, outc))

rm(outc, covs)

write.csv(data, file.path('data', 'DNAm_data.csv')) 

# ==== Plots ===================================================================
library(ggplot2)
library(pastaDaGg)

array_color_map <- c("450k" = "#FFC3CB", "EPICv1" = "#58aaa1", "EPICv2" = "#1c4b75")

dens_by_period <- ggplot(data, aes(x = cpg, colour = Period)) +
  geom_density() +
  labs(x = target_cpg, y = "Density", colour = "Period", 
       title = paste(target_cpg, "distribution by wave (Generation R)"))

dens_by_array <- ggplot(data, aes(x = cpg, color = Array, fill = Array)) + 
  geom_density(alpha = 0.3) +
  facet_grid(~ Period) +
  scale_color_manual(values = array_color_map) +
  scale_fill_manual(values = array_color_map) +
  labs(x = target_cpg, 
       y = "Density",
       title = paste(target_cpg, "distribution by wave and array (Generation R)"),
       color = "Array", fill = "Array") +
  theme_bw() +
  theme(legend.position = "bottom")

trajectories <- spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
          color = "Array", # split_by = "Array",
          title = paste(target_cpg, "trajectories (Generation R)"))

trajectories_smoking <- spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
                          color = "smoke", # split_by = "Array",
                          title = paste(target_cpg, "trajectories (Generation R)"))

trajectories_by_sex <- spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
                          color = "smoke", split_by = "sex",
                          title = paste(target_cpg, "trajectories (Generation R)"))

pdf(paste0(target_cpg,".pdf"), width = 12, height = 6)
dens_by_period
dens_by_array
trajectories
trajectories_smoking
trajectories_by_sex
dev.off()


