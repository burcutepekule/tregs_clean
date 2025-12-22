rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# ============================================================================
# READ PARAMETERS FROM CSV (WITH OPTIMIZED TREG PARAMETERS)
# ============================================================================

cat("Reading parameters with optimized Treg parameters...\n")
params_df = read.csv("./lhs_parameters_with_optimized_tregs.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

split_equal = function(x, n_chunks) {
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

# Get unique param_set_id and scenario_ind combinations
param_scenario_combos = params_df %>%
  select(param_set_id, scenario_ind) %>%
  distinct() %>%
  arrange(param_set_id, scenario_ind)

combo_ids = seq_len(nrow(param_scenario_combos)) - 1  # 0-indexed

chunks    = split_equal(combo_ids, n1)
loop_over = chunks[[n2]]

cat("Processing chunk", n2, "of", n1, "\n")
cat("Parameter-scenario combinations:", min(loop_over), "-", max(loop_over), "\n\n")

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/scratch/gpfs/CMETCALF/sim_abm_optimized'
dir.create(dir_name_data, showWarnings = FALSE)

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
t_max      = 5000
num_reps   = 10

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

k_in  = 0.044
x0_in = 50

cat("Simulation parameters:\n")
cat("  t_max:", t_max, "\n")
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(combo_idx in loop_over){
  scenario_elapsed_total = 0

  # Get param_set_id and scenario_ind for this combination
  param_set_id_use = param_scenario_combos$param_set_id[combo_idx + 1]  # +1 for 1-indexing
  scenario_ind_use = param_scenario_combos$scenario_ind[combo_idx + 1]

  # Get parameters for this combination
  param_set_use = params_df %>%
    filter(param_set_id == param_set_id_use, scenario_ind == scenario_ind_use)

  if (nrow(param_set_use) == 0) {
    cat("Warning: No parameters found for param_set_id", param_set_id_use,
        "scenario_ind", scenario_ind_use, "\n")
    next
  }

  # Reconstruct scenario from stored scenario_ind
  # Note: This requires reconstructing the scenarios_df from DLL_datagen_abm.R
  scenarios_df = expand.grid(
    control         = c(0),
    sterile         = c(0),
    allow_tregs     = c(0, 1),
    randomize_tregs = c(0),
    macspec_on      = c(0),
    ros_level       = c(0, 1),
    pat_level       = c(1, 2, 3)
  )
  scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
  scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==1))
  scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==0))

  scenarios_df_ctrl = expand.grid(
    control         = c(1),
    sterile         = c(0),
    allow_tregs     = c(0),
    randomize_tregs = c(0),
    macspec_on      = c(0),
    ros_level       = c(0),
    pat_level       = c(1)
  )
  scenarios_df=rbind(scenarios_df_ctrl, scenarios_df)

  sterile         = scenarios_df[scenario_ind_use,]$sterile
  allow_tregs     = scenarios_df[scenario_ind_use,]$allow_tregs
  randomize_tregs = scenarios_df[scenario_ind_use,]$randomize_tregs
  macspec_on      = scenarios_df[scenario_ind_use,]$macspec_on
  control         = scenarios_df[scenario_ind_use,]$control
  ros_level       = scenarios_df[scenario_ind_use,]$ros_level
  pat_level       = scenarios_df[scenario_ind_use,]$pat_level

  source("./MISC/ASSIGN_PARAMETERS.R")

  cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
             ' - scenario ', scenario_ind_use))

  scenario_start_time = Sys.time()

  longitudinal_df_keep = c()

  # ========================================================================
  # RUN SIMULATION WITH C++ ACCELERATION
  # ========================================================================
  source("./MISC/RUN_REPS_CPP_ABM_PAMPS.R")

  scenario_end_time = Sys.time()
  scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
  scenario_elapsed_total = scenario_elapsed_total + scenario_elapsed
  cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))

  saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                       '_control_',control,
                                       '_sterile_',sterile,
                                       '_macspec_',macspec_on,
                                       '_tregs_',allow_tregs,
                                       '_ros_level_',ros_level,
                                       '_pat_level_',pat_level,
                                       '_trnd_',randomize_tregs,
                                       '_optimized.rds'))

  cat(sprintf(' - %.1f seconds in total ✓\n', scenario_elapsed_total))
}
