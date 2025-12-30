rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(purrr)
library(stringr)

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
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

# args   = commandArgs(trailingOnly = TRUE)
# n1     = as.integer(args[1])
# n2     = as.integer(args[2])
# 
# chunks    = split_equal(0:max(params_df$param_set_id), n1)
# loop_over = chunks[[n2]]

loop_over = 60800
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
# params_df$recruitment_rate_danger = 0.01 # for 801

# cat("Processing chunk", n2, "of", n1, "\n")
cat("Parameter sets:", min(loop_over), "-", max(loop_over), "\n\n")

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/Users/burcutepekule/Desktop/sim_abm'
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

scenarios_df_exist = list.files(dir_name_data, pattern = "^longitudinal_df_param_set_id_.*\\.rds$") %>%
  str_match("sterile_([\\d.]+)_macspec_([\\d.]+)_tregs_([\\d.]+)_ros_level_([\\d.]+)_pat_level_([\\d.]+)_trnd_([\\d.]+)\\.rds") %>%
  as_tibble(.name_repair = "minimal") %>%
  set_names(c("full_match", "sterile", "macspec_on", "allow_tregs", "ros_level", "pat_level", "randomize_tregs")) %>%
  dplyr::select(-full_match) %>%
  distinct()

scenarios_df_exist = scenarios_df_exist %>%
  mutate(
    sterile = as.numeric(sterile),
    macspec_on = as.numeric(macspec_on),
    allow_tregs = as.numeric(allow_tregs),
    ros_level = as.numeric(ros_level),
    pat_level = as.numeric(pat_level),
    randomize_tregs = as.numeric(randomize_tregs)
  )

ros_levels_in = sort(unique(scenarios_df_exist$ros_level))

# scenarios_df_exist = expand.grid(
#   sterile         = c(0),
#   allow_tregs     = c(0),
#   randomize_tregs = c(0),
#   macspec_on      = c(0),
#   ros_level       = c(0,0.05,0.1,0.25,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5:10), # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
#   pat_level       = c(1:10)
# )

scenarios_df_all = expand.grid(
  sterile         = c(0),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0),
  macspec_on      = c(0),
  ros_level       = ros_levels_in, 
  pat_level       = c(1:15)
)

scenarios_df = anti_join(scenarios_df_all, scenarios_df_exist)
dim(scenarios_df)

# # DOESN'T MAKE SENSE TO RUN THIS
# scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
# scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==1))
# scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on>0 & allow_tregs == 1 & randomize_tregs==0))
# scenarios_df_ctrl = expand.grid(
#   sterile         = c(0),
#   allow_tregs     = c(0),
#   randomize_tregs = c(0),
#   macspec_on      = c(0),
#   ros_level       = c(0), # ros_level=0 makes this control
#   pat_level       = c(1)
# )
# scenarios_df=rbind(scenarios_df_ctrl, scenarios_df)

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over)*nrow(scenarios_df)*num_reps, "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  scenario_elapsed_total = 0
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
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

