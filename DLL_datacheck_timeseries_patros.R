rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

path         = "/Users/burcutepekule/Desktop/sim_abm/"
level_in     = 1
param_id_vec = readRDS(paste0('evo_selected_level_',level_in,'.rds'))
rep_ind_vec  = 0:9
alpha_plot   = 1/length(rep_ind_vec)

for(param_id in param_id_vec){
  
  results_0_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds'))
  results_0_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_2_trnd_0.rds'))
  results_0_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_3_trnd_0.rds'))
  results_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_1_trnd_0.rds'))
  results_1_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_2_trnd_0.rds'))
  results_1_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_3_trnd_0.rds'))
  results_2_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_1_trnd_0.rds'))
  results_2_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_2_trnd_0.rds'))
  results_2_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_3_trnd_0.rds'))
  
  results = rbind(
    results_0_1, results_0_2, results_0_3,
    results_1_1, results_1_2, results_1_3,
    results_2_1, results_2_2, results_2_3
  )
  
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

