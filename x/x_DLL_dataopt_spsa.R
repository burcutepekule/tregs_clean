rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
# library(cowplot)
# library(av)
# library(ggplot2)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

loop_over_all   = 1300
loop_over_all   = 60800

params_df       = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)

loop_over = loop_over_all
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

params_df_treg = params_df[param_names]
param_bounds = data.frame(
  lower = c(0.001, 0.001, 0.001, 0.750, 0.001),
  upper = c(0.120, 0.500, 0.250, 1.000, 1.000),
  row.names = param_names
)
params_df_treg = sapply(param_names, function(p) {
  runif(1, param_bounds[p, "lower"], param_bounds[p, "upper"])
})

params_df[param_names]=params_df_treg

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
num_reps   = 1
t_max      = 5000000
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
act_radius_PAMPs = 1

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

### replicate as you want to parallelize same scenario multiple times
scenarios_df = expand.grid(
  sterile         = c(0),
  allow_tregs     = c(1), # THIS NEEDS TO BE 1 AT ALL TIMES FOR TREG OPTIMIZATION!
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = c(3), # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  pat_level       = c(8)
)

# for (i in 1:4){
#   scenarios_df = rbind(scenarios_df, scenarios_df)
# }
# 
# dim(scenarios_df)[1]

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")

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
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1]) # this is not necessary right now 
n2     = as.integer(args[2]) # just use n2 to start another optimization

chunks = split_equal(1:nrow(scenarios_df), n1)
loop_over_sc = chunks[[n2]]


# ============================================================================
# DETECTION SETTINGS
# ============================================================================
collapse_threshold = 30
collapse_duration  = 25
collapse_rate      = 0.75

spsa_update_interval = 1000 # run this much to see any change by the parameter set

success_threshold_e  = 150*0.75
success_threshold_p  = 10
success_duration     = 50 # because max active_age_limit=30 and this ensures 2 cycles of the max value
success_rate         = 0.75
# success_threshold_e= 150*0.95 
# success_threshold_p= 10
# success_duration   = 250
# success_rate       = 0.95
# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  results = c()
  
  # for (scenario_ind in 1:dim(scenarios_df)[1]){
  # for (scenario_ind in 13){
  for (scenario_ind in loop_over_sc){
    
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
    # RUN SIMULATION WITH C++ ACCELERATION AND MACROPHAGE SPECIFICITY
    # ========================================================================
    source("./MISC/RUN_REPS_CPP_ABM_PAMPS_SPSA.R")
    
    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
    
    results = rbind(results, longitudinal_df_keep)
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}



