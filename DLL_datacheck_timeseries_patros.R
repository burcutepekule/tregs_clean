rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Define the ros and pat value ranges
ros_vals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
pat_vals = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

path         = "/Users/burcutepekule/Desktop/sim_abm/"
# level_in     = 1
# param_id_vec = readRDS(paste0('evo_selected_level_',level_in,'.rds'))]
param_id_vec = 0
rep_ind_vec  = 0:4
alpha_plot   = 1/length(rep_ind_vec)

for(param_id in param_id_vec){

  # Dynamically read all RDS files for this parameter set
  results_list = list()
  for (ros in ros_vals) {
    for (pat in pat_vals) {
      var_name = paste0("results_", ros, "_", pat)
      file_path = paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_', ros, '_pat_level_', pat, '_trnd_0.rds')
      results_list[[var_name]] = readRDS(file_path)
    }
  }

  # Combine all results
  results = do.call(rbind, results_list)
  
  results = results %>% dplyr::filter(rep_id %in% rep_ind_vec)

  
  # variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
  variables = c("epithelial_score")
  
  data_long = results %>%
    dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = alpha_plot, linewidth = 1) +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  
  ggsave(
    filename = paste0("./timeseries/level_",level_in,"_",param_id,"_",variables[1],".png"),
    plot = p,
    width = 14,
    height = 6,
    dpi = 300,
    bg='white'
  )
  
  variables = c('pathogen')
  
  data_long = results %>%
    dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = alpha_plot, linewidth = 1) +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() + scale_y_log10() +
    labs(title = "Pathogen Abundance", x = "Time", y = "Count", color = "Agent")
  
  # print(p)
  
  ggsave(
    filename = paste0("./timeseries/level_",level_in,"_",param_id,"_",variables[1],".png"),
    plot = p,
    width = 14,
    height = 6,
    dpi = 300,
    bg='white'
  )
  
}

