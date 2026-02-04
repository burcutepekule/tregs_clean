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
# params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
params_df = readRDS("./params_df.rds")
cat("Loaded", nrow(params_df), "parameter sets\n")

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
# CHUNKS
# ============================================================================

# PAY ATTENTION HERE!
source('./MISC/LOAD_PAT_LEVELS.R') # loads pat_level_vectors
# loop_over = as.numeric(names(pat_level_vectors))
loop_over = 43646
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
# ============================================================================
# HARDCODED PARAMETERS FOR PARAM_SET_ID 43646
# ============================================================================

cat("Using hardcoded parameters for param_set_id 43646...\n")

# params_df = data.frame(
#   param_set_id = 43646,
#   rate_leak_pathogen_injury = 0.5,
#   rate_leak_commensal_injury = 0.25,
#   rate_leak_commensal_baseline = 0.05,
#   epith_recovery_chance = 0.05,
#   th_ROS_microbe = 0.1088882,
#   th_ROS_epith_injury = 0.6158904,
#   diffusion_speed_DAMPs = 0.01648709,
#   diffusion_speed_PAMPs = 0.07892915,
#   diffusion_speed_SAMPs = 0.06675038,
#   diffusion_speed_ROS = 0.06533034,
#   add_ROS = 0.442411,
#   add_DAMPs = 0.09279615,
#   add_SAMPs = 0.06391013,
#   add_PAMPs = 0.3279585,
#   ros_decay = 0.1147807,
#   DAMPs_decay = 0.03812688,
#   SAMPs_decay = 0.1390501,
#   PAMPs_decay = 0.01297679,
#   activation_threshold_danger = 0.9478865,
#   activation_threshold_SAMPs = 0.4551507,
#   activity_engulf_M0_baseline = 0.05,
#   activity_engulf_M1_baseline = 0.2884559,
#   activity_ROS_M1_baseline = 0.4665152,
#   cc_phagocyte = 11,
#   active_age_limit = 7,
#   treg_discrimination_efficiency = 0.2402565,
#   recruitment_rate_danger = 0.2169343,
#   activity_engulf_M2_baseline = 0.2884559
# )
# cat("Loaded 1 parameter set (hardcoded)\n\n")

params_df$activity_engulf_M0_baseline = 0.00

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/scratch/gpfs/CMETCALF/sim_abm'
# dir_name_data = '/Users/burcutepekule/Desktop/sim_abm_local' # PAY ATTENTION HERE!
dir.create(dir_name_data, showWarnings = TRUE)

cat("Output directory:", dir_name_data, "\n\n")

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
source('./MISC/LOAD_FIXED_PARAMS.R') #num_reps = 10
# t_max      = 25
num_reps   = 100

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
# ===================== When running allow_tregs=1, optimized ================
scenarios_df = c()
for (param_id_in in loop_over){
  
  ### IF YOU WANNA INLCUDE OPT_IDX 0 CASE
  scenarios_df = rbind(scenarios_df, expand.grid(
    param_set_id    = param_id_in,
    sterile         = c(0),
    allow_tregs     = c(0), # PAY ATTENTION HERE!
    randomize_tregs = c(0),
    macspec_on      = c(0),
    ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
    pat_level       = pat_level_vectors[[as.character(param_id_in)]],
    overwrite       = c(1),
    diffusion_speed_SAMPs          = 0.1, # numbers so that it doesn't give NA or Inf somewhere
    add_SAMPs                      = 0.5, # numbers so that it doesn't give NA or Inf somewhere
    SAMPs_decay                    = 0.2, # numbers so that it doesn't give NA or Inf somewhere
    treg_discrimination_efficiency = 1, # numbers so that it doesn't give NA or Inf somewhere
    activation_threshold_SAMPs     = 0.25, # numbers so that it doesn't give NA or Inf somewhere
    opt_index                      = 0 # numbers so that it doesn't give NA or Inf somewhere
  ))
  
  # params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
  # # params_opt   = params_opt %>% dplyr::filter(pat_level==2) # for 630
  # params_opt   = params_opt[order(params_opt$mean_pct_above_threshold_min, params_opt$pat_level, decreasing = TRUE),]
  # params_opt   = na.omit(params_opt[c(1,3,6,14),])
  # print(params_opt)
  # for (ind_opt in 1:dim(params_opt)[1]){
  #   scenarios_df = rbind(scenarios_df, expand.grid(
  #     param_set_id    = param_id_in,
  #     sterile         = c(0),
  #     allow_tregs     = c(1), # PAY ATTENTION HERE!
  #     randomize_tregs = c(0),
  #     macspec_on      = c(0),
  #     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  #     pat_level       = pat_level_vectors[[as.character(param_id_in)]],
  #     overwrite       = c(1),
  #     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
  #     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
  #     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
  #     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
  #     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
  #     opt_index                      = ind_opt
  #   ))
  # }
}

dim(scenarios_df) 
cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

chunks       = split_equal(1:nrow(scenarios_df), n1)
loop_over_sc = chunks[[n2]]

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for (scenario_ind in loop_over_sc){
  
  scenario_elapsed_total = 0
  
  param_set_id_use = scenarios_df[scenario_ind,]$param_set_id
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
  sterile         = scenarios_df[scenario_ind,]$sterile
  allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
  randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
  macspec_on      = scenarios_df[scenario_ind,]$macspec_on
  ros_level       = scenarios_df[scenario_ind,]$ros_level
  pat_level       = scenarios_df[scenario_ind,]$pat_level
  overwrite_in    = scenarios_df[scenario_ind,]$overwrite
  opt_index       = scenarios_df[scenario_ind,]$opt_index
  
  param_names = c("diffusion_speed_SAMPs",
                  "add_SAMPs",
                  "SAMPs_decay",
                  "treg_discrimination_efficiency",
                  "activation_threshold_SAMPs")
  
  param_set_use[param_names] = c(scenarios_df[scenario_ind,]$diffusion_speed_SAMPs, 
                                 scenarios_df[scenario_ind,]$add_SAMPs, 
                                 scenarios_df[scenario_ind,]$SAMPs_decay, 
                                 scenarios_df[scenario_ind,]$treg_discrimination_efficiency,
                                 scenarios_df[scenario_ind,]$activation_threshold_SAMPs) #optimized for 68752
  
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
                                       '_overwrite_',overwrite_in,
                                       '_optidx_',opt_index,
                                       '.rds'))
  
}
cat(sprintf(' - %.1f seconds in total ✓\n', scenario_elapsed_total))

