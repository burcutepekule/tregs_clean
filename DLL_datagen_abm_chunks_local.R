rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
split_equal = function(x, n_chunks) {
  if (n_chunks == 1) {
    return(list(`1` = x))
  }
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

# ============================================================================
# FIRST HALF
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])
loop_over = split_equal(params_df$param_set_id, n1)[[n2]]

loop_over = 87503
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

# params_df[param_names] = c(0.027, 0.375, 0.016, 0.8, 0.221) #optimized for 68752
params_df[param_names] = c(0.009, 0.196, 0.415, 0.959, 0.958) #optimized for 87503
# params_df[param_names] = c(0.108, 0.107, 0.246, 0.905, 0.956) #optimized for 87503
# params_df[param_names] = c(0.059, 0.03, 0.241, 0.835, 0.132) #optimized for 25003
# params_df[param_names] = c(0.043, 0.238, 0.076, 0.944, 0.918) #optimized for 25003
# params_df[param_names] = c(0.054, 0.178, 0.380, 0.956, 0.111) #optimized for 25003

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/Users/burcutepekule/Desktop/sim_abm'
dir.create(dir_name_data, showWarnings = TRUE)

cat("Output directory:", dir_name_data, "\n\n")

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
source('./MISC/LOAD_FIXED_PARAMS.R')

colnames_insert = c('epithelial_score',
                    'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                    'commensal','pathogen','treg_resting','treg_active',
                    'C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')


cat("Simulation parameters:\n")
cat("  t_max:", t_max, "\n")
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n\n")

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================
scenarios_df = expand.grid(
  sterile         = c(0),
  allow_tregs     = c(0, 1), # PAY ATTENTION HERE! 
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  pat_level       = c(1, 2, 5, 7, 10, 12, 15, 20, 25),
  overwrite       = c(1)
  # ros_level       = seq(10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  # pat_level       = c(60,70,80,90,100)
)

dim(scenarios_df) 
cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

n3     = as.integer(args[3])
n4     = as.integer(args[4])

chunks        = split_equal(1:nrow(scenarios_df), n3)
loop_over_sc = chunks[[n4]]

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  scenario_elapsed_total = 0
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
  for (scenario_ind in loop_over_sc){
    
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    ros_level       = scenarios_df[scenario_ind,]$ros_level
    pat_level       = scenarios_df[scenario_ind,]$pat_level
    overwrite_in    = scenarios_df[scenario_ind,]$overwrite
    
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
    scenario_elapsed_total = scenario_elapsed_total + scenario_elapsed
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
    
    saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                         '_sterile_',sterile,
                                         '_macspec_',macspec_on,
                                         '_tregs_',allow_tregs,
                                         '_ros_level_',ros_level,
                                         '_pat_level_',pat_level,
                                         '_trnd_',randomize_tregs,
                                         '.rds'))
    
  }
  cat(sprintf(' - %.1f seconds in total ✓\n', scenario_elapsed_total))
}

