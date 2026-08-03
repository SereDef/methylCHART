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


# ==== Plots ===================================================================
library(ggplot2)
library(pastaDaGg)

target_pheno  <- c("mom_age_bin", "mom_age") # "smoke"
target_tag <- "mom_age"


array_color_map <- c("450k" = "#FFC3CB", "EPICv1" = "#58aaa1", "EPICv2" = "#1c4b75")

dens_by_period <- ggplot(testset, aes(x = Zcpg, colour = Period)) +
  geom_density() +
  labs(x = target_cpg, y = "Density", colour = "Period", 
       title = paste(target_cpg, "distribution by wave (Generation R)"))

dens_by_array <- ggplot(testset, aes(x = Zcpg, color = Array, fill = Array)) + 
  geom_density(alpha = 0.3) +
  facet_grid(~ Period) +
  scale_color_manual(values = array_color_map) +
  scale_fill_manual(values = array_color_map) +
  labs(x = target_cpg, 
       y = "Density",
       title = paste(target_cpg, "(Z) distribution by wave and array (Generation R)"),
       color = "Array", fill = "Array") +
  theme_bw() +
  theme(legend.position = "bottom")

trajectories <- spaghetti(testset,  x="Age", y = "Zcpg", id="IDC", interactive = FALSE, 
                          color = "Array", # split_by = "Array",
                          title = paste(target_cpg, "(Z) trajectories (Generation R)")) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")

pdf(paste0(target_tag, '_', target_cpg,"_z.pdf"), width = 12, height = 6)
dens_by_period
dens_by_array
trajectories
for (pheno in target_pheno) {
  traj_by <- spaghetti(testset,  x="Age", y = "Zcpg", id="IDC", interactive = FALSE, 
                       color = pheno, 
                       title = paste(target_cpg, " (Z) trajectories (Generation R) by", pheno))  +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")
  
  if (is.numeric(data[[pheno]])) {
    traj_by <- traj_by + scale_color_viridis_c(option = "magma")
  }
  
  # spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
  #           color = target_pheno, split_by = "sex",
  #           title = paste(target_cpg, "trajectories (Generation R)"))
  print(traj_by)
}
dev.off()



rm(dens_by_period, dens_by_array, trajectories, trajectories_smoking, trajectories_by_sex)
