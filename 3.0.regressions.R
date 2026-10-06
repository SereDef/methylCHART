# library(readr)
library(lmerTest)

# Read input data
cpg_data <- read.csv(file.path('data', 'DNAm_data.csv'), row.names = 1)
cov_data <- read.csv(file.path('data', 'cov_data.csv'), row.names = 1)

data = merge(cov_data, cpg_data, on='Sample_ID')

Z_test <- readr::read_csv("results/warpedBLR/results/Z_cpgdata_test.csv")
Z_train <- readr::read_csv("results/warpedBLR/results/Z_cpgdata_train.csv")

Z_data <- rbind(Z_test, Z_train) |>
  dplyr::select(-observations) |>
  dplyr::rename_with(~ifelse(.x == "subject_ids", "Sample_ID", paste0("Z", .x)))

alldata <- merge(data, Z_data, by='Sample_ID')

target_cpgs = jsonlite::read_json('cpg_targets.json', simplifyVector = FALSE, null = "null")

outp_regr <- list()

for (cpg in names(target_cpgs)) {
  if (length(target_cpgs[[cpg]]$pheno) > 0) {
    
    print(target_cpgs[[cpg]]$label)
    
    for (pheno in target_cpgs[[cpg]]$pheno) {
      
      if (!paste0('Z',cpg) %in% names(alldata)) {
        outp_regr[[cpg]][[pheno]]$raw <- NULL
        outp_regr[[cpg]][[pheno]]$norm <- NULL
        next
      }
      
      mod_form <- paste(cpg, '~', pheno, '+ Age + sex + (1|IDC)')
      
      fit1 <- lmer(as.formula(mod_form), data = alldata)
      fit2 <- lmer(as.formula(paste0('Z', mod_form)), data = alldata)
                   
       outp_regr[[cpg]][[pheno]]$raw <- capture.output(summary(fit1))
       outp_regr[[cpg]][[pheno]]$norm <- capture.output(summary(fit2))
    }
  }
}

library(jsonlite)

write_json(outp_regr, "model_summary.json", pretty = TRUE, auto_unbox = TRUE)

testset$Zcpg_bined <- cut(
  testset$Zcpg,
  breaks = c(-Inf, -2, -1, 1, 2, Inf),
  labels = c("< -2", "-2 to -1", "-1 to 1", "1 to 2", "> 2"),
  right = TRUE   # (-Inf,-2], (-2,-1], ...
)


quick_stats <- function(f, d = testset) {
  print(aggregate(as.formula(f), data=d, mean))
}

for (pheno in target_pheno[1]) {
  
  quick_stats(paste("cpg ~", pheno))
  quick_stats(paste("Zcpg ~", pheno))
  
  quick_stats(paste("cpg ~", pheno, '+ Array'))
  quick_stats(paste("Zcpg ~", pheno, '+ Array'))
  
  quick_stats(paste("cpg ~", pheno, '+ Period'))
  quick_stats(paste("Zcpg ~", pheno, '+ Period'))
}

aggregate(I(mom_age_bin == "> 35") ~ Zcpg_bined, data=testset, mean)
aggregate(I(mom_age_bin == "> 35") ~ Zcpg_bined + Period, data=testset, mean)



slow_stats <- function(rhs, d = testset) {
  fit1 <- lmer(as.formula(paste('cpg ~', rhs)), data = d)
  fit2 <- lmer(as.formula(paste('Zcpg ~', rhs)), data = d)
  
  print(summary(fit1))
  print(summary(fit2))
}

slow_stats('mom_age_bin + Age + sex + (1|IDC)')
slow_stats('mom_age + Age + sex + (1|IDC)')
