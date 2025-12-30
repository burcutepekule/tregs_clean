rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(deSolve)
library(ggplot2)
library(cowplot)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/ODE_SYSTEM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
loop_over = 6001
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
# params_df$activity_engulf_M1_baseline = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2',
                    'epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                    'commensal','pathogen','treg_resting','treg_active',
                    'C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
plot_on    = 0
plot_every = 0
t_max      = 2000

grid_size       = 25
n_phagocytes    = round(grid_size*grid_size*0.20)
n_tregs         = round(grid_size*grid_size*0.20)
n_commensals_lp = 20
max_total_phagocytes = round(grid_size*grid_size*0.80)

injury_percentage = 60

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1
max_cell_value_PAMPs = 1

act_radius_ROS   = 1
act_radius_treg  = 1
act_radius_DAMPs = 1
act_radius_SAMPs = 1

# Logistic function parameters (for epithelial injury calculation)
k_in  = 0.044
x0_in = 50

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

scenarios_df = expand.grid(
  sterile         = c(0),
  allow_tregs     = c(0),
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = c(0,1,2,3,4,5,6,7,8,9,10),
  pat_level       = c(1,2,3,4,5,6,7,8,9,10)
)

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  scenario_elapsed_total = 0
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  results = c()
  
  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    ros_level       = scenarios_df[scenario_ind,]$ros_level
    pat_level       = scenarios_df[scenario_ind,]$pat_level

    source("./MISC/ASSIGN_PARAMETERS.R")

    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))

    # Track timing for this scenario
    scenario_start_time = Sys.time()

    longitudinal_df_keep = c()

    # ========================================================================
    # RUN ODE SIMULATION
    # ========================================================================
    source("./MISC/RUN_REPS_ODE_PAMPS.R")

    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
    scenario_elapsed_total = scenario_elapsed_total + scenario_elapsed
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))

    results = rbind(results, longitudinal_df_keep)

  }
  cat(sprintf('  Total for param set %d: %.1f seconds ✓\n', param_set_id_use, scenario_elapsed_total))
  
  variables = c("epithelial_score")
  
  data_long = results %>%
    dplyr::select(t, tregs_on, ros_level, pat_level, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value") 
  
  p_e = ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
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
    dplyr::select(t, tregs_on, ros_level, pat_level, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p_p = ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
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
    filename = paste0("./more_ros/param_set_id_",param_set_id_use,".png"),
    plot = p_all,
    width = 22,
    height = 14,
    dpi = 300,
    bg='white'
  )
  
  ggsave(
    filename = paste0("./more_ros/param_set_id_",param_set_id_use,"_e.png"),
    plot = p_e,
    width = 14,
    height = 12,
    dpi = 300,
    bg='white'
  )
  
  ggsave(
    filename = paste0("./more_ros/param_set_id_",param_set_id_use,"_p.png"),
    plot = p_p,
    width = 14,
    height = 12,
    dpi = 300,
    bg='white'
  )
}

cat("\nODE data generation complete!\n")
