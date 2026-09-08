# library(readr)

# data <- readr::read_csv(file.path('data', 'DNAm_data.csv')) 

target_cpg = "cg00528572"


Z_test <- readr::read_csv(paste0("results/", target_cpg, "_warpedBLR/results/Z_cpgdata_test.csv"))

testset <- data[Z_test$observations+1, ]
testset$Zcpg <- Z_test$cpg

rm(Z_test)

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

library(lmerTest)

slow_stats <- function(rhs, d = testset) {
  fit1 <- lmer(as.formula(paste('cpg ~', rhs)), data = d)
  fit2 <- lmer(as.formula(paste('Zcpg ~', rhs)), data = d)
  
  print(summary(fit1))
  print(summary(fit2))
}

slow_stats('mom_age_bin + Age + sex + (1|IDC)')
slow_stats('mom_age + Age + sex + (1|IDC)')
