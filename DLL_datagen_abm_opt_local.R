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
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

loop_over     = c(50002)
params_df     = params_df %>% dplyr::filter(param_set_id %in% loop_over)
df_opt_rnd_95 = readRDS(paste0('./df_opt_rnd_95_',loop_over,'_use.rds'))

# chunks = split_equal(1:dim(df_opt_rnd_95)[1], n1)
# loop_over_param_inds = chunks[[n2]]
loop_over_param_inds = 1:dim(df_opt_rnd_95)[1]

for (i_opt in loop_over_param_inds){
  
  print(paste0('At i_opt', i_opt))
  params_df[param_names] = df_opt_rnd_95[i_opt,]
  title_opt = paste(df_opt_rnd_95[i_opt,], collapse = "_")
  
  cat("Processing chunk", n2, "of", n1, "\n")
  cat("Parameter sets:", min(loop_over), "-", max(loop_over), "\n\n")
  
  # ============================================================================
  # SETUP OUTPUT DIRECTORY
  # ============================================================================
  
  dir_name_data <<- '/Users/burcutepekule/Desktop/sim_abm_local'  # note the <<- for global assignment
  dir.create(dir_name_data, showWarnings = TRUE)
  
  cat("Output directory:", dir_name_data, "\n\n")
  
  colnames_insert = c('epithelial_score',
                      'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                      'commensal','pathogen','treg_resting','treg_active',
                      'C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')
  
  source('./MISC/LOAD_FIXED_PARAMS.R')
  num_reps   = 3
  
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
    allow_tregs     = c(1), # PAY ATTENTION HERE! 
    randomize_tregs = c(0),
    macspec_on      = c(0),
    # ros_level       = c(0,1,3,5,10), # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
    ros_level       = c(1,3,5,10), # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
    pat_level       = c(1,2,5,10,100)
  )
  
  cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
  cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")
  
  # ============================================================================
  # MAIN SIMULATION LOOP
  # ============================================================================
  chunks        = split_equal(1:nrow(scenarios_df), n1)
  loop_over_sc = chunks[[n2]]
  
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
                                           '_trnd_',randomize_tregs,'.rds'))
      
    }
    cat(sprintf(' - %.1f seconds in total ✓\n', scenario_elapsed_total))
  }
  # source('~/Dropbox/tregs_clean/DLL_datacheck_timeseries_patros.R')
}
