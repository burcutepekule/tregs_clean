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
# CHUNKS
# ============================================================================

# PAY ATTENTION HERE!
source('./MISC/LOAD_PAT_LEVELS.R') # loads pat_level_vectors
# loop_over = as.numeric(names(pat_level_vectors))
loop_over = c(17, 23, 55, 58, 88)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

params_df$activity_engulf_M0_baseline = 0.00
params_df$activity_engulf_M1_baseline = 0.05
params_df$activity_engulf_M2_baseline = 0.05
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
  params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
  params_opt   = params_opt[params_opt$mean_pct_above_threshold_e>=.9,]
  params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
  params_opt   = na.omit(params_opt[1:2,])
  print(params_opt)
  for (ind_opt in 1:dim(params_opt)[1]){
    scenarios_df = rbind(scenarios_df, expand.grid(
      param_set_id    = param_id_in,
      sterile         = c(0),
      allow_tregs     = c(1), # PAY ATTENTION HERE!
      randomize_tregs = c(0),
      macspec_on      = c(0),
      ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
      pat_level       = pat_level_vectors[[as.character(param_id_in)]],
      overwrite       = c(0,1),
      diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
      add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
      SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
      treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
      activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
      opt_index                      = ind_opt
    ))
  }
}

dim(scenarios_df) 
# param_names = c("diffusion_speed_SAMPs",
#                 "add_SAMPs",
#                 "SAMPs_decay",
#                 "treg_discrimination_efficiency",
#                 "activation_threshold_SAMPs")
# param_id_in  = 68752
# params_opt   = read.table(paste0("./merged_",param_id_in,".txt"),
#                              header = FALSE,
#                              col.names = c("param_set_id" ,"rep", "pat_level", "ros_level",
#                                            param_names, "min_e", "min_p","mean_e","mean_p",
#                                            "pct_above_threshold_e","pct_above_threshold_p"))
# 
# params_opt = params_opt %>% dplyr::filter(pct_above_threshold_e==1 & pct_above_threshold_p==1)
# params_opt = distinct(params_opt[c("param_set_id","diffusion_speed_SAMPs","add_SAMPs","SAMPs_decay","treg_discrimination_efficiency","activation_threshold_SAMPs")])
# 
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       =  c(1, 5, 7, 10, 15, 20, 25, 30),
#     overwrite       = c(0, 1),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }

# param_id_in  = 92000
# params_opt   = readRDS(paste0('./summary_df_045_',param_id_in,'_use.rds'))
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = c(1, 5, 7,seq(8, 12, 1)),
#     overwrite       = c(0),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }

# param_id_in  = 30000
# params_opt   = readRDS(paste0('./summary_df_050_',param_id_in,'_use.rds'))
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = c(1, seq(2, 5, 0.5)),
#     overwrite       = c(0),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }

# for(param_id_in in c(250, 47500, 66250, 78750, 97250)){
#   params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
#   params_opt   = params_opt[params_opt$mean_pct_above_threshold_e==1,]
#   params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
#   params_opt   = params_opt[1:3,]
#   params_opt   = na.omit(params_opt) # in case less than 10
#   for (ind_opt in 1:dim(params_opt)[1]){
#     scenarios_df = rbind(scenarios_df, expand.grid(
#       param_set_id    = param_id_in,
#       sterile         = c(0),
#       allow_tregs     = c(1), # PAY ATTENTION HERE!
#       randomize_tregs = c(0),
#       macspec_on      = c(0),
#       ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#       pat_level       = pat_level_vectors[[as.character(param_id_in)]],
#       overwrite       = c(0,1),
#       diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#       add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#       SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#       treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#       activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#       opt_index                      = ind_opt
#     ))
#   }
# }
# 
# param_id_in  = 67250
# params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
# params_opt   = params_opt[params_opt$mean_pct_above_threshold_e>=0.4,]
# params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
# params_opt   = params_opt[1,]
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = pat_level_vectors[[as.character(param_id_in)]],
#     overwrite       = c(0,1),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }
# 
# param_id_in  = 80750
# params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
# params_opt   = params_opt[params_opt$mean_pct_above_threshold_e>=0.4,]
# params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
# params_opt   = params_opt[1,]
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = pat_level_vectors[[as.character(param_id_in)]],
#     overwrite       = c(0,1),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }
# 
# param_id_in  = 90250
# params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
# params_opt   = params_opt[params_opt$mean_pct_above_threshold_e>=0.4,]
# params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
# params_opt   = params_opt[1,]
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = pat_level_vectors[[as.character(param_id_in)]],
#     overwrite       = c(0,1),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }

# param_id_in  = 73750
# params_opt   = readRDS(paste0('./summary_df_10rep_',param_id_in,'_use.rds'))
# params_opt   = params_opt[params_opt$mean_pct_above_threshold_e>=0.4,]
# params_opt   = params_opt[order(params_opt$pat_level, params_opt$ros_level, decreasing = TRUE),]
# params_opt   = params_opt[1:3,]
# for (ind_opt in 1:dim(params_opt)[1]){
#   scenarios_df = rbind(scenarios_df, expand.grid(
#     param_set_id    = param_id_in,
#     sterile         = c(0),
#     allow_tregs     = c(1), # PAY ATTENTION HERE!
#     randomize_tregs = c(0),
#     macspec_on      = c(0),
#     ros_level       = seq(0,10,1), # MAX 10! 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#     pat_level       = pat_level_vectors[[as.character(param_id_in)]],
#     overwrite       = c(0,1),
#     diffusion_speed_SAMPs          = params_opt[ind_opt,]$diffusion_speed_SAMPs,
#     add_SAMPs                      = params_opt[ind_opt,]$add_SAMPs,
#     SAMPs_decay                    = params_opt[ind_opt,]$SAMPs_decay,
#     treg_discrimination_efficiency = params_opt[ind_opt,]$treg_discrimination_efficiency,
#     activation_threshold_SAMPs     = params_opt[ind_opt,]$activation_threshold_SAMPs,
#     opt_index                      = ind_opt
#   ))
# }

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

