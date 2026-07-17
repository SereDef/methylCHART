# library(readr)

# data <- readr::read_csv(file.path('data', 'DNAm_data.csv')) 

Z_test <- readr::read_csv("BLR_warped/results/Z_cpgdata_test.csv")

testset <- data[Z_test$observations+1, ]
testset$Zcpg <- Z_test$cpg

rm(Z_test)

aggregate(cpg ~ smoke, data = testset, mean)
aggregate(Zcpg ~ smoke, data = testset, mean)

aggregate(cpg ~ smoke + Array, data = testset, mean)
aggregate(Zcpg ~ smoke + Array, data = testset, mean)

aggregate(cpg ~ smoke + Period, data = testset, mean)
aggregate(Zcpg ~ smoke + Period, data = testset, mean)

# ==== Plots ===================================================================
library(ggplot2)
library(pastaDaGg)

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
                          title = paste(target_cpg, "(Z) trajectories (Generation R)"))

trajectories_smoking <- spaghetti(testset,  x="Age", y = "Zcpg", id="IDC", interactive = FALSE, 
                                  color = "smoke", # split_by = "Array",
                                  title = paste(target_cpg, "(Z) trajectories (Generation R)"))

trajectories_by_sex <- spaghetti(testset,  x="Age", y = "Zcpg", id="IDC", interactive = FALSE, 
                                 color = "smoke", split_by = "sex",
                                 title = paste(target_cpg, "(Z) trajectories (Generation R)"))

pdf(paste0(target_cpg,"_z.pdf"), width = 12, height = 6)
dens_by_period
dens_by_array
trajectories
trajectories_smoking
trajectories_by_sex
dev.off()

rm(dens_by_period, dens_by_array, trajectories, trajectories_smoking, trajectories_by_sex)
