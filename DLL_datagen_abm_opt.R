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

loop_over = c(50002)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
opt_df    = readRDS(paste0('./df_opt_rnd_95_',loop_over,'_use.rds'))

param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/scratch/gpfs/CMETCALF/sim_abm'
dir.create(dir_name_data, showWarnings = TRUE)

cat("Output directory:", dir_name_data, "\n\n")

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
num_reps   = 5

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
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n\n")

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================
scenarios_df = expand.grid(
  sterile         = c(0),
  allow_tregs     = c(1), # TO CHECK OPTIMIZATION ALWAYS 1 
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = c(0, 1, seq(2, 12, 0.5)), # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  pat_level       = c(1, 2, seq(3, 10, 0.5))
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
    
    for(optidx in 1:dim(opt_df)[1]){# == OPTIMIZED FOR 60800
      cat("Processing chunk", optidx, "of", dim(opt_df)[1], "\n")

      params_insert = opt_df[optidx,]
      param_set_use[param_names] = params_insert ## OPTIMIZED FOR 60800_A - WORKS BEAUTIFULLY FOR _A!
      title_opt = paste(params_insert, collapse = "_")
      
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
                                           '_optidx_',optidx,
                                           '.rds'))
    }
  }
  cat(sprintf(' - %.1f seconds in total ✓\n', scenario_elapsed_total))
}

