rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(cowplot)
library(av)
library(ggplot2)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# reslevel_in     = 0
# loop_over_all   = readRDS(paste0('evo_selected_reslevel_',reslevel_in,'.rds'))
# loop_over_all   = 1401
# loop_over_all   = readRDS('evo_selected_js_based.rds')

loop_over_all   = 8100
params_df_keep  = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
params_df_96101 = params_df_keep %>% dplyr::filter(param_set_id==96101)
param_names = c("param_set_id",
                "diffusion_speed_SAMPs", 
                "add_SAMPs", 
                "SAMPs_decay", 
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")
params_df_96101 = params_df_96101[param_names]
params_df_96101_long = params_df_96101 %>% 
  pivot_longer(cols = -param_set_id, 
               names_to = "parameter", 
               values_to = "value")

# args   = commandArgs(trailingOnly = TRUE)
# n1     = as.integer(args[1])
# n2     = as.integer(args[2])
# 
# chunks    = split_equal(loop_over_all, n1)
# loop_over = chunks[[n2]]

loop_over = loop_over_all

### OPTIMIZED VIA SPSPA

opts_df = data.frame(
  try = 1:41,
  diffusion_speed_SAMPs  = c(0.010, 0.010, 0.010, 0.010, 0.001, 0.014,
                             0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.0096, 0.001, 0.001, 0.002, 0.010, 0.001, 0.017, 0.001, 0.0029, 0.001,
                             0.001, 0.001, 0.001, 0.001, 0.0047, 0.001, 0.0034, 0.001, 0.0262, 0.018, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
                             0.001),
  add_SAMPs              = c(0.111, 0.500, 0.327, 0.289, 0.294, 0.102,
                             0.500, 0.001, 0.0973, 0.1214, 0.001, 0.139, 0.165, 0.4356, 0.2187, 0.1325, 0.1707, 0.1375, 0.3915, 0.2346, 0.271, 0.0474, 0.1469, 0.2924,
                             0.500, 0.184, 0.1916, 0.4206, 0.1359, 0.4985, 0.500, 0.4507, 0.500, 0.1657, 0.1251, 0.001, 0.0904, 0.079, 0.479, 0.4267,
                             0.001),
  SAMPs_decay            = c(0.111, 0.340, 0.298, 0.221, 0.222, 0.009,
                             0.043, 0.001, 0.0332, 0.0083, 0.001, 0.0713, 0.0326, 0.1718, 0.0312, 0.0169, 0.1292, 0.1645, 0.0575, 0.0592, 0.0084, 0.0162, 0.0103, 0.3284,
                             0.0798, 0.0243, 0.2031, 0.0646, 0.0184, 0.0743, 0.193, 0.1657, 0.0741, 0.0414, 0.103, 0.001, 0.0434, 0.0125, 0.022, 0.0453,
                             0.001),
  treg_discrimination_efficiency = c(0.932, 0.245, 0.785, 1.000, 0.406, 0.746,
                                     0.2626, 0.2126, 0.0026, 0.1348, 0.9682, 0.5246, 0.6645, 0.9396, 0.5877, 0.705, 0.1995, 0.1622, 0.3327, 0.8118, 0.304, 0.5444, 0.271, 0.2941,
                                     0.6543, 1.000, 0.7031, 0.703, 0.3545, 0.3331, 0.2568, 0.2184, 0.5787, 0.918, 0.3025, 0.3256, 0.497, 0.4521, 0.5528, 0.511,
                                     0.001),
  activation_threshold_SAMPs     = c(0.091, 0.278, 0.241, 0.010, 0.385, 0.324,
                                     0.8898, 0.2143, 0.8438, 0.5258, 0.1793, 0.4777, 0.7109, 0.5281, 0.6685, 0.5577, 0.2972, 0.3045, 0.7611, 0.4656, 0.9732, 0.7348, 0.8202, 0.4898,
                                     0.6776, 0.4607, 0.4189, 0.4793, 0.4029, 0.4403, 0.3733, 0.4786, 0.8989, 0.8388, 0.4006, 0.2585, 0.5691, 0.8107, 0.5909, 0.6946,
                                     0.2389))
# Check for command line argument
args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  try_indices = as.integer(args[1])
} else {
  try_indices = 38:41
}

for (try_in in try_indices){
  
  params_df = params_df_keep %>% dplyr::filter(param_set_id %in% loop_over)
  
  opts_df_use = opts_df %>% dplyr::filter(try==try_in)
  
  params_df$diffusion_speed_SAMPs          = opts_df_use$diffusion_speed_SAMPs
  params_df$add_SAMPs                      = opts_df_use$add_SAMPs
  params_df$SAMPs_decay                    = opts_df_use$SAMPs_decay
  params_df$treg_discrimination_efficiency = opts_df_use$treg_discrimination_efficiency
  params_df$activation_threshold_SAMPs     = opts_df_use$activation_threshold_SAMPs
  
  # ============================================================================
  # SETUP OUTPUT DIRECTORY
  # ============================================================================
  
  colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2',
                      'epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                      'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                      'commensal','pathogen','treg_resting','treg_active',
                      'C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')
  
  # ============================================================================
  # FIXED PARAMETERS (not in CSV)
  # ============================================================================
  num_reps   = 5
  t_max      = 5000
  
  plot_on    = 0
  plot_every = 10
  if(plot_on==1){
    dir_name_data = '/Users/burcutepekule/Desktop/gif_out'
    dir.create(dir_name_data, showWarnings = FALSE)
    cat("Output directory:", dir_name_data, "\n\n")
    dir_name_frames = paste0(dir_name_data,'/frames')
    dir.create(dir_name_frames, showWarnings = FALSE)
  }
  grid_size       = 25
  n_phagocytes    = round(grid_size*grid_size*0.20)
  n_tregs         = round(grid_size*grid_size*0.20)
  n_commensals_lp = 20
  max_total_phagocytes = round(grid_size*grid_size*0.80)
  
  injury_percentage = 60
  max_level_injury  = 5
  
  max_cell_value_ROS   = 1
  max_cell_value_DAMPs = 1
  max_cell_value_SAMPs = 1
  max_cell_value_PAMPs = 1
  
  ## PLOTTING
  lim_ROS  = max_cell_value_ROS
  lim_DAMP = max_cell_value_DAMPs
  lim_SAMP = max_cell_value_SAMPs
  lim_PAMP = max_cell_value_PAMPs
  ## PLOTTING
  
  act_radius_ROS   = 1
  act_radius_treg  = 1
  act_radius_DAMPs = 1
  act_radius_SAMPs = 1
  
  # Logistic function parameters (for epithelial injury calculation)
  k_in  = 0.044
  x0_in = 50
  
  cat("Simulation parameters:\n")
  cat("  t_max:", t_max, "\n")
  cat("  num_reps:", num_reps, "\n")
  cat("  grid_size:", grid_size, "x", grid_size, "\n")
  cat("  n_phagocytes:", n_phagocytes, "\n")
  cat("  n_tregs:", n_tregs, "\n\n")
  
  # ============================================================================
  # SCENARIO DEFINITIONS
  # ============================================================================
  
  scenarios_df = expand.grid(
    # sterile         = c(0, 1),
    sterile         = c(0),
    allow_tregs     = c(0, 1),
    # randomize_tregs = c(0, 1),
    randomize_tregs = c(0),
    # macspec_on      = c(0, 1, 2)
    macspec_on      = c(0),
    ros_level       = c(1, 2), 
    pat_level       = c(1, 2, 3)
  )
  # DOESN'T MAKE SENSE TO RUN THIS
  scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
  scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==1))
  scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==0))
  scenarios_df_ctrl = expand.grid(
    # sterile         = c(0, 1),
    sterile         = c(0),
    allow_tregs     = c(0),
    randomize_tregs = c(0),
    macspec_on      = c(0),
    ros_level       = c(0),
    pat_level       = c(1)
  )
  scenarios_df=rbind(scenarios_df_ctrl, scenarios_df)
  
  cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
  cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")
  
  # ============================================================================
  # MAIN SIMULATION LOOP
  # ============================================================================
  
  # # ===== to check
  # scenarios_df = scenarios_df %>% dplyr::filter(ros_level==2 & allow_tregs==1)
  # # ===== to check
  
  for(param_set_id_use in loop_over){
    param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
    results = c()
    
    for (scenario_ind in 1:dim(scenarios_df)[1]){
      
      sterile         = scenarios_df[scenario_ind,]$sterile
      allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
      randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
      macspec_on      = scenarios_df[scenario_ind,]$macspec_on
      ros_level       = scenarios_df[scenario_ind,]$ros_level
      pat_level       = scenarios_df[scenario_ind,]$pat_level
      
      source("./MISC/ASSIGN_PARAMETERS.R")
      
      cat(paste0('[', Sys.time(), '] Processing optimal try ',try_in,' param set ', param_set_id_use,
                 ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))
      
      # Track timing for this scenario
      scenario_start_time = Sys.time()
      
      longitudinal_df_keep = c()
      
      # ========================================================================
      # RUN SIMULATION WITH C++ ACCELERATION AND MACROPHAGE SPECIFICITY
      # ========================================================================
      source("./MISC/RUN_REPS_CPP_ABM_PAMPS.R")
      
      scenario_end_time = Sys.time()
      scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
      
      results = rbind(results, longitudinal_df_keep)
      
      cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
    }
    
    variables = c("epithelial_score")
    
    data_long = results %>%
      dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value") 
    
    p_e = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
      geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
      facet_grid(tregs_on ~ ros_level + pat_level, labeller = label_both) +
      scale_color_manual(values = agent_colors) +
      theme_minimal() +
      theme(
        # legend.position = "none", #turn legend off
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
      )+
      labs(title = "", x = "Time", y = "Count", color = "Agent") 
    # labs(title = paste0(variables, " dynamics"), x = "Time", y = "Count", color = "Agent") 
    
    variables = c("pathogen")
    
    data_long = results %>%
      dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    p_p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
      geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
      facet_grid(tregs_on ~ ros_level + pat_level, labeller = label_both) +
      scale_color_manual(values = agent_colors) +
      theme_minimal() +
      theme(
        # legend.position = "none", #turn legend off
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 16),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
      )+
      labs(title = "", x = "Time", y = "Count", color = "Agent") 
    # labs(title = paste0(variables, " dynamics"), x = "Time", y = "Count", color = "Agent") 
    
    graphics.off()
    title_grob = ggdraw() + draw_label(paste0("param_set_id: ",param_set_id_use), fontface = "bold", size = 20)
    p_all = plot_grid(title_grob, p_e, p_p, ncol = 1, rel_heights = c(0.05, 1, 1))
    
    ggsave(
      filename = paste0("./more_ros/param_set_id_",param_set_id_use,"_try_",try_in,".png"),
      plot = p_all,
      width = 22,
      height = 14,
      dpi = 300,
      bg='white'
    )
    
  }
}
