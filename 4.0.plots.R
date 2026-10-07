# ==== Plots ===================================================================
library(ggplot2)
library(pastaDaGg)

model_names <- list.files('results')

for (model_name in model_names) {
  
  model_desc <- strsplit(model_name, "_")[[1]] 
  names(model_desc) <- c('algorithm', 'likelihood', 'sample', 'age_model')
  
  # Read input data
  cpg_data <- read.csv(file.path('data', 'dnam_data.csv'), row.names = 1)
  cov_data <- read.csv(file.path('data', 'pheno_data.csv'), row.names = 1)
  
  data <- merge(cov_data, cpg_data, on='Sample_ID')
  
  model_path <- file.path('results', model_name)
  
  Z_test <- readr::read_csv(file.path(model_path, 'results', 'Z_test.csv'))
  Z_train <- readr::read_csv(file.path(model_path, 'results', 'Z_train.csv'))
  
  Z_data <- rbind(Z_test, Z_train) |>
    dplyr::select(-observations) |>
    dplyr::rename_with(~ifelse(.x == "subject_ids", "Sample_ID", paste0("Z", .x)))
  
  # Merge, keep only the rows used in the model train/test (e.g. GENR only)
  alldata <- merge(data, Z_data, by='Sample_ID', all.x = FALSE, all.y = TRUE) |>
    dplyr:: mutate(Period = factor(Period, 
      levels = c("Birth", "Age3.5", "Age5", "Age7", "Age9", "Age13", "Age15", "Age17", "Age24")))
  
  target_cpgs = jsonlite::read_json('cpg_targets.json', simplifyVector = FALSE, null = "null")
  
  plots_location <- file.path(model_path, 'plots')
  
  make_plots <- function(cpg, d = alldata, z = FALSE, outp = plots_location, 
                         age_var = model_desc[['age_model']],
                         array_color_map = c("450k" = "#FFC3CB", "EPICv1" = "#58aaa1", "EPICv2" = "#1c4b75"), 
                         cohort_color_map = c('GENR' = '#2c3c94', 'ALSPAC'='#9D021f')
                         # split_color_map <- c('train' = '', 'test'='')
                         ) {
    label <- target_cpgs[[cpg]]$label
    pheno <- target_cpgs[[cpg]]$pheno
    
    z_label <- if (z) ' Z' else ''
    z_lab <- gsub(' ', '', z_label)
    cpg_var <- if (z) paste0('Z',cpg) else cpg
    
    if (!cpg_var %in% names(d)) {
      cli::cli_inform('{cpg_var} [{label}] is missing!')
      return(NULL)
    }
    
    dens_by_period <- ggplot(d, aes(x = .data[[cpg_var]], colour = Period)) +
      geom_density() +
      labs(x = cpg, y = "Density", colour = "Period", 
           title = paste0(cpg, z_label,' [',label,'] distribution by wave'))
    
    ggsave(file.path(outp, paste0('dens_', cpg, z_lab, '_by_period.png')), 
           dens_by_period, width = 8, height = 4, units = "in", dpi = 300)
    
    dens_by_array <- ggplot(d, aes(x = .data[[cpg_var]], color = Array, fill = Array)) + 
      geom_density(alpha = 0.3) +
      facet_grid(~ Period) +
      scale_color_manual(values = array_color_map) +
      scale_fill_manual(values = array_color_map) +
      labs(x = cpg, 
           y = "Density",
           title = paste0(cpg, z_label,' [',label,'] distribution by wave and array'),
           color = "Array", fill = "Array") +
      theme_bw() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(outp, paste0('dens_', cpg, z_lab, '_by_array.png')), 
           dens_by_array, width = 12, height = 4, units = "in", dpi = 300)
    
    dens_by_cohort <- ggplot(d, aes(x = .data[[cpg_var]], color = Cohort, fill = Cohort)) + 
      geom_density(alpha = 0.3) +
      facet_grid(~ Array) +
      scale_color_manual(values = cohort_color_map) +
      scale_fill_manual(values = cohort_color_map) +
      labs(x = cpg, 
           y = "Density",
           title = paste0(cpg, z_label,' [',label,'] distribution by array and cohort'),
           color = "Cohort", fill = "Cohort") +
      theme_bw() +
      theme(legend.position = "bottom")
    
    ggsave(file.path(outp, paste0('dens_', cpg, z_lab, '_by_cohort.png')), 
           dens_by_cohort, width = 8, height = 4, units = "in", dpi = 300)
    
    # --- Trajectories -----------------------------------------------------------
    traj_colors <- c('Array', 'Cohort', unlist(pheno))
    
    for (p in traj_colors) {
      traj_by <- spaghetti(d,  x = age_var, y = cpg_var, id="IDC", interactive = FALSE, 
                           color = p, 
                           title = paste0(cpg, z_label,' [',label,'] trajectories')) +
        theme(legend.position = "bottom")
      
      if (p == 'Array') {
        traj_by <- traj_by + scale_color_manual(values = array_color_map)
      } else if (p == 'Cohort') {
        traj_by <- traj_by + scale_color_manual(values = cohort_color_map)
      } else if (is.numeric(d[[p]])) {
        traj_by <- traj_by + scale_color_viridis_c(option = "magma", na.value = NA)
      } else {
        traj_by <- traj_by + scale_color_discrete(na.translate = FALSE)
      }
      
      if (z) {
        traj_by <- traj_by +
          geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black")
      }
      
      ggsave(file.path(outp, paste0('traj_', cpg, z_lab, '_by_', p,'.png')), traj_by,
             width = 9, height = 7, units = "in", dpi = 300)
    }
  
  }
  
  for (cpg in names(cpg_data)[-1]) {
    print(cpg)
    make_plots(cpg, z = FALSE)
    make_plots(cpg, z = TRUE)
  }

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


