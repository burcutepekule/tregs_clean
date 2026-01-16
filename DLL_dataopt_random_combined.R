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
loop_over = c(62500, 67250, 73750, 80750)
# ============================================================================
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

params_df_treg = params_df[param_names]
param_bounds = data.frame(
  lower = c(0.001, 0.001, 0.001, 0.750, 0.001),
  upper = c(0.120, 0.500, 0.500, 1.000, 1.000),
  row.names = param_names
)

# Don't set initial parameters - random search will sample them
# Just keep param_bounds for the search

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

scenarios_df_1 = expand.grid(
  param_set_id    = c(62500),
  sterile         = c(0),
  allow_tregs     = c(1), # THIS NEEDS TO BE 1 AT ALL TIMES!
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = c(5:10),
  pat_level       = c(8.5,9,9.5,10),
  overwrite       = c(1)
)
# scenarios_df_2 = expand.grid(
#   param_set_id    = c(73750, 80750),
#   sterile         = c(0),
#   allow_tregs     = c(1), # THIS NEEDS TO BE 1 AT ALL TIMES!
#   randomize_tregs = c(0),
#   macspec_on      = c(0),
#   ros_level       = c(3:10),
#   pat_level       = c(9,10),
#   overwrite       = c(1)
# )
# scenarios_df_3 = expand.grid(
#   param_set_id    = c(67250),
#   sterile         = c(0),
#   allow_tregs     = c(1), # PAY ATTENTION HERE! 
#   randomize_tregs = c(0),
#   macspec_on      = c(0),
#   ros_level       = c(4:10),
#   pat_level       = c(3, 3.5, 4, 4.5, 5),
#   overwrite       = c(1)
# )
# scenarios_df = rbind(scenarios_df_1, scenarios_df_2, scenarios_df_3)
scenarios_df = rbind(scenarios_df_1)
dim(scenarios_df) 
rownames(scenarios_df)=1:dim(scenarios_df)[1]
scenarios_df = scenarios_df[rep(seq_len(nrow(scenarios_df)), each = 17), ]
dim(scenarios_df) 

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

chunks        = split_equal(1:nrow(scenarios_df), n1)
loop_over_sc = chunks[[n2]]

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/scratch/gpfs/CMETCALF/sim_opt_random'
# dir_name_data = '/Users/burcutepekule/Desktop/sim_opt_random'
dir.create(dir_name_data, showWarnings = TRUE)

cat("Output directory:", dir_name_data, "\n\n")

# ============================================================================
# DETECTION SETTINGS
# ============================================================================
success_threshold_e = 5
success_threshold_p = 10
success_duration    = 150
success_rate        = 0.95
max_iterations      = 10000 # Number of random samples to try
# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for (scenario_ind in loop_over_sc){
  results = c()

  param_set_id_use = scenarios_df[scenario_ind,]$param_set_id
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
  sterile         = scenarios_df[scenario_ind,]$sterile
  allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
  randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
  macspec_on      = scenarios_df[scenario_ind,]$macspec_on
  ros_level       = scenarios_df[scenario_ind,]$ros_level
  pat_level       = scenarios_df[scenario_ind,]$pat_level
  overwrite_in    = scenarios_df[scenario_ind,]$overwrite
  
  source("./MISC/ASSIGN_PARAMETERS.R")
  
  cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
             ' ros ', ros_level, ' pat ', pat_level, 
             ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))
  
  # Track timing for this scenario
  scenario_start_time = Sys.time()
  
  # ========================================================================
  # RUN RANDOM SEARCH OPTIMIZATION
  # ========================================================================
  source("./MISC/RUN_REPS_CPP_ABM_PAMPS_RANDOM.R")
  
  scenario_end_time = Sys.time()
  scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
  
  cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  
}

