rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(cowplot)
library(av)
library(ggplot2)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

params_df    = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
loop_over    = c(38002)
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
num_reps   = 10
t_max      = 2000
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
n_phagocytes    = round(grid_size * grid_size * 0.20)
n_tregs         = round(grid_size * grid_size * 0.20)
n_commensals_lp = 20

injury_percentage = 60
max_level_injury  = 5

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

lim_ROS  = max_cell_value_ROS
lim_DAMP = max_cell_value_DAMPs
lim_SAMP = max_cell_value_SAMPs

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
  sterile         = c(0, 1),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0, 1),
  macspec_on      = c(0, 1)
)
# DOESN'T MAKE SENSE TO RUN THIS
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on ==1 & allow_tregs == 1 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on ==1 & allow_tregs == 1 & randomize_tregs==0))
scenarios_df_ctrl = expand.grid(
  control         = c(1),
  sterile         = c(0, 1),
  allow_tregs     = c(0),
  randomize_tregs = c(0),
  macspec_on      = c(0)
)
scenarios_df=rbind(scenarios_df_ctrl, scenarios_df)

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over) * nrow(scenarios_df) * num_reps, "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
results = c()
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  # for (scenario_ind in 1:nrow(scenarios_df)){
  for (scenario_ind in 1){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    control         = scenarios_df[scenario_ind,]$control
    
    source("./MISC/ASSIGN_PARAMETERS.R")
    
    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))

    # Track timing for this scenario
    scenario_start_time = Sys.time()

    longitudinal_df_keep = c()

    # ========================================================================
    # RUN SIMULATION WITH C++ ACCELERATION AND MACROPHAGE SPECIFICITY
    # ========================================================================
    source("./MISC/RUN_REPS_CPP_ABM.R")

    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
    
    results = rbind(results, longitudinal_df_keep)
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}

variables = c("epithelial_score")

data_long = results %>%
  dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = .1, linewidth = 1) +
  facet_grid(randomize_tregs ~ control + macspec_on + sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")

print(p)
