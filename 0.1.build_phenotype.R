use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.5'
.libPaths(use_library)

library(dplyr)
# library(haven)

# === INFO =====================================================================
# > GENR `DNAm_selection_file.sav` stored in the data subfolder is the selection 
#   file for Release 4 (got this from Isabel)
# > ALSPAC phenotype, link and samplesheet files are taken from the TEMPO project
#   updated on 19 Sept 2027 (B3361), the data subfolder 

# === Generation R =============================================================

genr <- haven::read_sav(file.path('data', 
                                  'DNAm_selection_file.sav')) |>
  # Some clean-up
  mutate(
    Cohort = 'GENR',
    Period = factor(Period, levels = c("Birth", "Age5", "Age9", "Age13", "Age17")),
    Array = as.factor(Batch),
    Batch = interaction(Period, Array, drop = TRUE), 
    Age2 = ifelse(Period == "Birth", (Gestational_Age - 40) / 52.18, Age))

summary(genr)

# Other covariates / phenotypes of interest ------------------------------------

covs <- haven::read_sav(file.path('data',
                                  'Family_GeneralData_Pregnancy_20251111.sav')) |>
  dplyr::transmute(IDM,
                   Sex = haven::as_factor(GENDERPREG), # 1 = boy 2 = girl
                   mom_age = AGE_M_Birth, # years
                   mom_age_bin = as.factor(ifelse(mom_age > 35, '> 35', '=< 35'))
  ) 

msmk <- haven::read_sav(file.path('data',
                                  'Mother_Smoking_Pregnancy_20260216.sav')) |>
  dplyr::transmute(IDM, mom_smoking = haven::as_factor(SMOKE_ALL))
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
  dplyr::transmute(IDC, wheez = as.factor(wh4y))

genr <- Reduce(function(x, y) merge(x, y, by = "IDM", all.x = TRUE), 
               list(genr, covs, msmk))
genr <- Reduce(function(x, y) merge(x, y, by = "IDC", all.x = TRUE), 
               list(genr, adhd, asma)) |> 
  mutate(IDC = as.character(IDC))

rm(covs, msmk, adhd, asma)

summary(genr)


# === ALSPAC ===================================================================

# Go fetch gestational age at birth from phenotype
alsp_pheno <- haven::read_sav('~/Alspac_Tempo/phenotype/B3361_Cecil_19Sept27.sav') |> 
  dplyr::transmute(pheno_id = paste0(cidB3361, qlet), 
                   Gestational_Age = as.numeric(bestgest)) |> 
  # Sex = factor(kz021, levels = c(1, 2), labels = c("boy", "girl"),
  # mom_age = mz028b,
  # mom_smoking = b032) |> 
  dplyr::left_join( 
    haven::read_sav('~/Alspac_Tempo/link_file/OmicsIDs_B3361_20Oct23.sav') |>  
      dplyr::transmute(pheno_id = paste0(cidB3361, qlet),
                       dnam_id = paste0(dnam_epic450_g0_g1, qlet)),
    by = 'pheno_id')

alsp <- read.csv('~/Alspac_Tempo/methylation/data/samplesheet/samplesheet.csv') |>
  dplyr::filter(qlet %in% c('A', 'B')) |> # Exclude parents
  dplyr::transmute(
    Cohort = 'ALSPAC',
    dnam_id = paste0(dnam_epic450_g0_g1, qlet),
    Sample_ID = Sample_Name,
    Array = factor(chip, levels = c("450k", "epic"), labels = c("450k", "EPICv1")),
    Period = factor(time_code, 
                    levels = c('cord', 'c43m',  'c61m', 'F7',  'F9',  'TF3',  'F17',  'F24'), 
                    labels = c('Birth','Age3.5','Age5', 'Age7','Age9','Age15','Age17','Age24')),
    Age = ifelse(Period =='Birth', 0, age),
    Sex = factor(Sex, levels = c("M", "F"), labels = c("boy", "girl"))
  ) |>
  dplyr::left_join(alsp_pheno, by = 'dnam_id') |> 
  dplyr::mutate(
    IDC = pheno_id, 
    Age2 = ifelse(Period == "Birth", (Gestational_Age - 40) / 52.18, Age)) |>
  # TMP mean impute missing age 
  group_by(Period) |>
  mutate(Age = ifelse(is.na(Age), median(Age, na.rm = TRUE), Age),
         Age2 = ifelse(is.na(Age2), median(Age2, na.rm = TRUE), Age2)) |>
  ungroup()

rm(alsp_pheno)

# === BIND TOGETHER ============================================================

summary(genr)
summary(alsp)

aggregate(Age2 ~ Period, summary, data = alsp)
aggregate(Age2 ~ Period, summary, data = genr)

data <- bind_rows(
  alsp |> select(Sample_ID, IDC, Cohort, Period, Array, Age, Age2, Sex, Gestational_Age),
  genr |> select(Sample_ID, IDC, Cohort, Period, Array, Age, Age2, Sex, Gestational_Age,
                 # Isolation_Lab, Run, Bridge
                 mom_age, mom_age_bin, mom_smoking, # the genr-only columns
                 cbcl_att, cbcl_adhd, cbcl_int, cbcl_ext, cbcl_tot, wheez)   
) |> 
  mutate(Period = factor(Period, 
                         levels = c("Birth", "Age3.5", "Age5", "Age7", "Age9", "Age13", "Age15", "Age17", "Age24")))

aggregate(Age2 ~ Period, summary, data = data)

aggregate(IDC ~ Cohort, data = data, FUN = dplyr::n_distinct)
aggregate(IDC ~ Period, data = data, FUN = dplyr::n_distinct)

table(data$Period, data$Array)

write.csv(data, file.path('data', 'pheno_data.csv'))
