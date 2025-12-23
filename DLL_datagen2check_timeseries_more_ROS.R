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

reslevel_in     = 0
loop_over_all   = readRDS(paste0('evo_selected_reslevel_',reslevel_in,'.rds'))
loop_over_all   = 96101

params_df       = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)

# df_results_keep = readRDS(paste0('./data_cpp_read_abm','_ros_vs_ctrl_patros','.rds'))
# df_results_pick = df_results_keep %>% dplyr::filter(param_set_id %in% loop_over_all)

# params_df_pick = params_df %>% dplyr::filter(param_set_id==2716)
# params_long_pick = params_df_pick %>%
#   pivot_longer(
#     cols = -param_set_id,  # keep param_set_id as identifier
#     names_to = "parameter",
#     values_to = "value"
#   )

# split_equal = function(x, n_chunks) {
#   split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
# }
# 
# args   = commandArgs(trailingOnly = TRUE)
# n1     = as.integer(args[1])
# n2     = as.integer(args[2])
# 
# chunks    = split_equal(loop_over_all, n1)
# loop_over = chunks[[n2]]

loop_over = loop_over_all
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

params_df$activity_engulf_M1_baseline = 0.2
params_df$activity_engulf_M2_baseline = 0.2

# ### make it equal to 96101 in terms of tregs? 
# params_df$activation_threshold_SAMPs = 0.0224
# params_df$add_SAMPs = 0.0243

# # ==== I mean they are not optimal of course, one can always find the optimal config for tregs to be most useful
# params_df$SAMPs_decay
# params_df$activation_threshold_SAMPs

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
  control         = c(0),
  # sterile         = c(0, 1),
  sterile         = c(0),
  allow_tregs     = c(0, 1),
  # randomize_tregs = c(0, 1),
  randomize_tregs = c(0),
  # macspec_on      = c(0, 1, 2)
  macspec_on      = c(0),
  ros_level       = c(1, 2), # for proper labelling
  pat_level       = c(1, 2, 3)
)
# DOESN'T MAKE SENSE TO RUN THIS
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==0))
scenarios_df_ctrl = expand.grid(
  control         = c(1),
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
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  results = c()
  
  for (scenario_ind in 1:13){
    
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    control         = scenarios_df[scenario_ind,]$control
    ros_level       = scenarios_df[scenario_ind,]$ros_level
    pat_level       = scenarios_df[scenario_ind,]$pat_level
    
    source("./MISC/ASSIGN_PARAMETERS.R")
    
    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
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
    dplyr::select(t, control, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value") %>%
    mutate(control = factor(control, levels = c(1, 0))) # to plot control=1 first
  
  p_e = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
    facet_grid(tregs_on ~ control + ros_level + pat_level, labeller = label_both) +
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
    dplyr::select(t, control, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")%>%
    mutate(control = factor(control, levels = c(1, 0))) # to plot control=1 first
  
  p_p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
    facet_grid(tregs_on ~ control + ros_level + pat_level, labeller = label_both) +
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
    filename = paste0("./more_ros/param_set_id_",param_set_id_use,"_reslevel_",reslevel_in,".png"),
    plot = p_all,
    width = 22,
    height = 14,
    dpi = 300,
    bg='white'
  )
  
}

