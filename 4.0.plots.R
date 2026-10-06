# ==== Plots ===================================================================
library(ggplot2)
library(pastaDaGg)

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

array_color_map <- c("450k" = "#FFC3CB", "EPICv1" = "#58aaa1", "EPICv2" = "#1c4b75")

plots_location <- './results/plots'

make_plots <- function(cpg, d = alldata, z = FALSE, outp = plots_location, 
                       array_color_map = c("450k" = "#FFC3CB", "EPICv1" = "#58aaa1", "EPICv2" = "#1c4b75")) {
  label <- target_cpgs[[cpg]]$label
  pheno <- target_cpgs[[cpg]]$pheno
  
  label2 <- if (z) ' Z' else ''
  cpg_var <- if (z) paste0('Z',cpg) else cpg
  
  if (!cpg_var %in% names(d)) {
    cli::cli_inform('{cpg_var} [{label}] is missing!')
    return(NULL)
  }
  
  dens_by_period <- ggplot(d, aes(x = .data[[cpg_var]], colour = Period)) +
    geom_density() +
    labs(x = cpg, y = "Density", colour = "Period", 
         title = paste0(cpg, label2,' [',label,'] distribution by wave (Generation R)'))
  
  ggsave(file.path(outp, paste0(cpg, label2, '_dens_by_period.png')), 
         dens_by_period)
  
  dens_by_array <- ggplot(d, aes(x = .data[[cpg_var]], color = Array, fill = Array)) + 
    geom_density(alpha = 0.3) +
    facet_grid(~ Period) +
    scale_color_manual(values = array_color_map) +
    scale_fill_manual(values = array_color_map) +
    labs(x = cpg, 
         y = "Density",
         title = paste0(cpg, label2,' [',label,'] distribution by wave and array (Generation R)'),
         color = "Array", fill = "Array") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(outp, paste0(cpg, label2, '_dens_by_array.png')), 
         dens_by_array)
  
  trajectories <- spaghetti(d,  x="Age", y = cpg_var, id="IDC", interactive = FALSE, 
                            color = "Array", # split_by = "Array",
                            title = paste0(cpg, label2,' [',label,'] trajectories (Generation R)'))
  
  if (z) {
    trajectories <- trajectories +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")
  }
  
  ggsave(file.path(outp, paste0(cpg, label2, '_trajectories.png')), 
         trajectories)
  
  if (length(pheno) > 0) {
    for (p in pheno) {
      traj_by <- spaghetti(d,  x="Age", y = cpg_var, id="IDC", interactive = FALSE, 
                           color = p, 
                           title = paste0(cpg, label2,' [',label,'] trajectories (Generation R)'))
      if (is.numeric(d[[p]])) {
        traj_by <- traj_by + scale_color_viridis_c(option = "magma")
      }
      
      if (z) {
        traj_by <- traj_by +
          geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")
      }
      ggsave(file.path(outp, paste0(cpg, label2, '_trajectories_', p,'.png')), traj_by)
    }
  } 

}

for (cpg in names(cpg_data)[-1]) {
  make_plots(cpg, z = TRUE)
}

# traj_by <- spaghetti(alldata,  x="Age", y = 'Zcg00528572', id="IDC", interactive = FALSE, 
#                      color = 'mom_age_bin',
#                      title = 'cg00528572 trajectories (Generation R)') +
#   geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")
# 
# saveRDS(traj_by, 'plot_for_alex.rds')
# pdf(paste0(target_tag, '_', target_cpg,".pdf"), width = 12, height = 6)
# dens_by_period
# dens_by_array
# trajectories
# for (pheno in target_pheno) {
#   traj_by <- spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
#                        color = pheno, 
#                        title = paste(target_cpg, "trajectories (Generation R) by", pheno)) 
#   
#   if (is.numeric(data[[pheno]])) {
#     traj_by <- traj_by + scale_color_viridis_c(option = "magma")
#   }
#   
#   # spaghetti(data,  x="Age", y = "cpg", id="IDC", interactive = FALSE, 
#   #           color = target_pheno, split_by = "sex",
#   #           title = paste(target_cpg, "trajectories (Generation R)"))
#   print(traj_by)
# }
# dev.off()


