rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

path         = "/Users/burcutepekule/Desktop/sim_abm/"
param_id_vec = 38002
rep_ind_vec  = 0:9
alpha_plot   = 1/length(rep_ind_vec)

control_pick         = c(1)
sterile_pick         = c(0, 1)
tregs_on_pick        = c(0)
macspec_on_pick      = c(0)
randomize_tregs_pick = c(0)

for(param_id in param_id_vec){
  # Control
  results_1_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
  results_1_1_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_1_sterile_1_macspec_0_tregs_0_trnd_0.rds'))
  
  # === STERILE
  # Test (with ROS)
  results_0_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
  results_0_0_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_0_macspec_0_tregs_1_trnd_0.rds'))
  results_0_0_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_0_macspec_0_tregs_1_trnd_1.rds'))
  # Test (with ROS), perfect macrophage
  results_0_0_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_0_macspec_1_tregs_0_trnd_0.rds'))
  
  # === PATHOGENIC
  # Test (with ROS)
  results_0_1_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_1_macspec_0_tregs_0_trnd_0.rds'))
  results_0_1_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_1_macspec_0_tregs_1_trnd_0.rds'))
  results_0_1_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_1_macspec_0_tregs_1_trnd_1.rds'))
  # Test (with ROS), perfect macrophage
  results_0_1_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_control_0_sterile_1_macspec_1_tregs_0_trnd_0.rds'))
  
  results = rbind(
    results_1_0_0_0_0, results_1_1_0_0_0,
    results_0_0_0_0_0, results_0_0_0_1_0,
    results_0_0_0_1_1, results_0_0_1_0_0,
    results_0_1_0_0_0, results_0_1_0_1_0,
    results_0_1_0_1_1, results_0_1_1_0_0
  )
  
  results = results %>% dplyr::filter(rep_id %in% rep_ind_vec)
  
  results = results %>% dplyr::filter(control %in% control_pick
                                      & sterile %in% sterile_pick 
                                      & tregs_on %in% tregs_on_pick
                                      & macspec_on %in% macspec_on_pick
                                      & randomize_tregs %in% randomize_tregs_pick)
  
  # variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
  variables = c("epithelial_score")
  
  data_long = results %>%
    dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = alpha_plot, linewidth = 1) +
    facet_grid(randomize_tregs ~ control + macspec_on + sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  
  # print(p)
  
  ggsave(
    filename = paste0("./",variables,"_",param_id,".png"),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300,
    bg='white'
  )
  
  variables = c('commensal')
  
  data_long = results %>%
    dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = alpha_plot, linewidth = 1) +
    facet_grid(randomize_tregs ~ control + macspec_on + sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  
  # print(p)
  
  ggsave(
    filename = paste0("./",variables,"_",param_id,".png"),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300,
    bg='white'
  )
  
  variables = c('pathogen')
  
  data_long = results %>%
    dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = alpha_plot, linewidth = 1) +
    facet_grid(randomize_tregs ~ control + macspec_on + sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  
  # print(p)
  
  ggsave(
    filename = paste0("./",variables,"_",param_id,".png"),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300,
    bg='white'
  )
}

