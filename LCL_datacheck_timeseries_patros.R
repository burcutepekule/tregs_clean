rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Define the ros and pat value ranges
# ros_vals = seq(0,5,0.25)
ros_vals = seq(0,5,0.25)
# ros_vals = 0
# pat_vals = seq(12,20,2)
# pat_vals = c(50,60,70,80,90,100,110,120,130,140,150)
pat_vals = c(0.1,2:49,50,60,70,80,90,100,110,120,130,140,150)

level_in     = 0 
path         = "/Users/burcutepekule/Desktop/sim_abm/"
param_id_vec = 801
rep_ind_vec  = 0:4
alpha_plot   = 1/length(rep_ind_vec)

for(param_id in param_id_vec){

  # # Check if output files already exist
  # file1 = paste0("/Users/burcutepekule/Desktop/timeseries_local/level_", level_in, "_", param_id, "_epithelial_score.png")
  # file2 = paste0("/Users/burcutepekule/Desktop/timeseries_local/level_", level_in, "_", param_id, "_pathogen.png")
  # 
  # if (file.exists(file1) && file.exists(file2)) {
  #   message("Skipping param_id ", param_id, " - files already exist")
  #   next
  # }
  
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
    filename = paste0("/Users/burcutepekule/Desktop/timeseries_local/level_",level_in,"_",param_id,"_",variables[1],".png"),
    plot = p,
    width = 24,
    height = 20,
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
    filename = paste0("/Users/burcutepekule/Desktop/timeseries_local/level_",level_in,"_",param_id,"_",variables[1],".png"),
    plot = p,
    width = 24,
    height = 20,
    dpi = 300,
    bg='white'
  )
  
}

