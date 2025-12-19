rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(cowplot)
library(av)
library(ggplot2)
library(changepoint)

# source('./DLL_datazanalyse_abm_agl_and_cycling.R')

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

params_df      = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
df_agl_cycling_osc  = readRDS('df_agl_cycling_osc.rds')
df_agl_cycling_nosc = readRDS('df_agl_cycling_nosc.rds')

# make two columns of ss_start out of one
df_agl_cycling_osc = df_agl_cycling_osc %>%
  dplyr::select(param_set_id, score_type, mean_ss_start) %>%
  pivot_wider(
    names_from = score_type,
    values_from = mean_ss_start,
    names_prefix = "mean_ss_start_"
  ) %>%
  rename(
    mean_ss_start_e = mean_ss_start_epithelium,
    mean_ss_start_p = mean_ss_start_pathogen
  ) %>%
  left_join(
    df_agl_cycling_osc %>% 
      dplyr::select(-score_type, -mean_ss_start) %>% 
      distinct(),
    by = "param_set_id"
  )

df_agl_cycling_nosc = df_agl_cycling_nosc %>%
  dplyr::select(param_set_id, score_type, mean_ss_start) %>%
  pivot_wider(
    names_from = score_type,
    values_from = mean_ss_start,
    names_prefix = "mean_ss_start_"
  ) %>%
  rename(
    mean_ss_start_e = mean_ss_start_epithelium,
    mean_ss_start_p = mean_ss_start_pathogen
  ) %>%
  left_join(
    df_agl_cycling_nosc %>% 
      dplyr::select(-score_type, -mean_ss_start) %>% 
      distinct(),
    by = "param_set_id"
  )

loop_over_df        = na.omit(rbind(df_agl_cycling_osc, df_agl_cycling_nosc))
loop_over_df        = loop_over_df %>% dplyr::filter(mean_ss_start_e<250 & mean_ss_start_p<250)  

loop_over_df_osc    = loop_over_df %>% dplyr::filter(oscillator_type=='oscillating')
loop_over_df_nosc   = loop_over_df %>% dplyr::filter(oscillator_type=='not oscillating')

n_check             = dim(loop_over_df_osc)[1]
loop_over_df        = na.omit(rbind(loop_over_df_osc[1:n_check,], loop_over_df_nosc[1:n_check,]))

loop_over_all        = loop_over_df$param_set_id
loop_over_label_type = loop_over_df$oscillator_type
loop_over_label_type = gsub('not oscillating', 'not_oscillating', loop_over_label_type)
loop_over_label_src  = loop_over_df$oscillating_source
loop_over_ss_e       = loop_over_df$mean_ss_start_e
loop_over_ss_p       = loop_over_df$mean_ss_start_p

## =============================================================================
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

loop_over    = loop_over_all
loop_over    = 3409
## =============================================================================

params_df    = params_df %>% dplyr::filter(param_set_id %in% loop_over)

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
t_max      = 750

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
  # sterile         = c(0),
  allow_tregs     = c(0, 1),
  # randomize_tregs = c(0),
  # macspec_on      = c(0, 1, 2),
  ros_level       = c(0, 1)
)
# DOESN'T MAKE SENSE TO RUN THIS
# scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
# scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==1))
# scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==0))
scenarios_df_ctrl = expand.grid(
  control         = c(1),
  # sterile         = c(0),
  allow_tregs     = c(0),
  # randomize_tregs = c(0),
  # macspec_on      = c(0),
  ros_level       = c(0)
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
  
  for (scenario_ind in c(2)){
  # for (scenario_ind in 1:5){
    
    sterile         = 0
    randomize_tregs = 0
    macspec_on      = 0
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    control         = scenarios_df[scenario_ind,]$control
    ros_level       = scenarios_df[scenario_ind,]$ros_level
    
    if(ros_level==1){
      param_set_use$add_ROS = 1
    }
    
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
    longitudinal_df_keep$ros_level = ros_level
    
    
    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
    
    results = rbind(results, longitudinal_df_keep)
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
  
  variables = c("epithelial_score")
  
  data_long_e = results %>%
    dplyr::select(t, control, tregs_on, ros_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p_e = ggplot(data_long_e, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
    facet_grid(ros_level ~ control + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = paste0(variables, " Dynamics"), x = "Time", y = "Count", color = "Agent")
  
  variables = c("pathogen")
  
  data_long_p = results %>%
    dplyr::select(t, control, tregs_on, ros_level, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p_p = ggplot(data_long_p, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = max(1/num_reps,.1), linewidth = 1) +
    facet_grid(ros_level ~ control + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = paste0(variables, " Dynamics"), x = "Time", y = "Count", color = "Agent")
  
  graphics.off()
  
  ind=which(loop_over==param_set_id_use)
  
  ss_plot_on_p = loop_over_ss_p[ind]
  ss_plot_on_e = loop_over_ss_p[ind]
  
  p_p = p_p + geom_vline(xintercept = ss_plot_on_p, linetype = "dashed", color = "magenta", linewidth=2)
  p_e = p_e + geom_vline(xintercept = ss_plot_on_e, linetype = "dashed", color = "magenta", linewidth=2)
  
  title_grob = ggdraw() + draw_label(paste0("param_set_id: ",param_set_id_use," agl: ",param_set_use$active_age_limit," type: ",loop_over_label_type[ind]," src: ",loop_over_label_src[ind]), fontface = "bold")
  p_all = plot_grid(title_grob, p_e, p_p, ncol = 1, rel_heights = c(0.1, 1, 1))

  ggsave(
    filename = paste0("./agl_cycling/",loop_over_label_type[ind],"_param_set_id_",param_set_id_use,".png"),
    plot = p_all,
    width = 26,
    height = 14,
    dpi = 300,
    bg='white'
  )
  
}

