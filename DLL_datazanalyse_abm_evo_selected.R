rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)
library(ggsignif)

jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = 25*25*0.05
M1_M2_diff     = 0
filter_control = 0
labels_on      = 1
score_type     = 'pathogen'
data_suffix    = '_ros_vs_ctrl_nopatlevel' #
data_suffix    = '_ros_vs_ctrl_patros' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
if(data_suffix == '_ros_vs_ctrl_nopatlevel'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else if(data_suffix == '_ros_vs_ctrl_patros'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==260) %>% dplyr::pull(Var1) # 260 = 10 reps per scenario, 13 scenarios x 2 times recording for epithelial and pathogen scores 
}
df_results = df_results_keep %>% filter(param_set_id %in% keep_param_id)
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500 # used to be 4500, just for simulation purposes to save time
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
num_params = length(unique(df_results$param_set_id))

df_comparisons = df_results %>% dplyr::select(
  param_set_id, injury_type, 
  starts_with("d_ros"), starts_with("mean_ros"))

df_comparisons_keep = distinct(df_comparisons)

# ============= HELPER FUNCTION: Define "infection under control" ====================================================
# Infection is under control when:
# 1. Epithelial health is approximately 150 (threshold: > 149)
# 2. Pathogen load is approximately 0 (threshold: < 1)
is_under_control = function(epithelial_health, pathogen_load) {
  return(epithelial_health > 149 & pathogen_load < 1)
}

# ============= STEP 1: Filter for evolutionary need for ROS ====================================================
# Rule 1: Without ROS (ros_0, treg0_), even the lowest pathogen level (pat1_) should:
#   a) Cause significant injury to epithelium (ros0_pat1 NOT under control)
#   b) Have higher pathogen burden than ros1_pat1 case (ros1_pat1 should be better)

df_step1 = df_comparisons_keep %>%
  dplyr::filter(injury_type == 'pathogenic') %>%
  dplyr::mutate(
    # Check if ros0_pat1_treg0 is NOT under control (infection causes damage)
    ros0_pat1_controlled = is_under_control(mean_ros0_pat1_treg0_e, mean_ros0_pat1_treg0_p),
    ros1_pat1_controlled = is_under_control(mean_ros1_pat1_treg0_e, mean_ros1_pat1_treg0_p),
    ros1_sigdiff_ros0    = (d_ros0_pat1_treg0_vs_ros1_pat1_treg0_p >= jsd_th & d_ros0_pat1_treg0_vs_ros1_pat1_treg0_e >= jsd_th )
  ) %>%
  dplyr::filter(!ros0_pat1_controlled & ros1_pat1_controlled)

cat("After Step 1 (ROS needed for pat1):", nrow(df_step1), "parameter sets\n")
evo_selected_reslevel_0_ids = sort(df_step1$param_set_id)
saveRDS(evo_selected_reslevel_0_ids, 'evo_selected_reslevel_0.rds')


df_step1_treg = df_step1[c('param_set_id',
                           'd_ros1_pat1_treg0_vs_ros1_pat1_treg1_e',
                           'd_ros2_pat1_treg0_vs_ros2_pat1_treg1_e',
                           'd_ros1_pat2_treg0_vs_ros1_pat2_treg1_e',
                           'd_ros2_pat2_treg0_vs_ros2_pat2_treg1_e',
                           'd_ros1_pat3_treg0_vs_ros1_pat3_treg1_e',
                           'd_ros2_pat3_treg0_vs_ros2_pat3_treg1_e')]

df_step1_96101   = df_step1 %>% dplyr::filter(param_set_id==96101) %>% dplyr::select(-injury_type)
params_long_pick = df_step1_96101 %>%
  pivot_longer(
    cols = -param_set_id,  # keep param_set_id as identifier
    names_to = "parameter",
    values_to = "value"
  )


# ============= STEP 2: Split into resistance levels ====================================================
# Identify two groups based on how much ROS is needed (still no tregs, treg0_):

# reslevel_1: ros1_ controls pat1_ but NOT pat2_
df_reslevel_1 = df_step1 %>%
  dplyr::mutate(
    ros1_pat2_worse = 
      # (mean_ros1_pat2_treg0_p>mean_ros1_pat1_treg0_p & d_ros1_pat1_treg0_vs_ros1_pat2_treg0_p >= jsd_th &
      # mean_ros1_pat2_treg0_e<mean_ros1_pat1_treg0_e & d_ros1_pat1_treg0_vs_ros1_pat2_treg0_e >= jsd_th)
      (mean_ros1_pat2_treg0_p>mean_ros1_pat1_treg0_p &
         mean_ros1_pat2_treg0_e<mean_ros1_pat1_treg0_e)
  ) %>%
  dplyr::filter(ros1_pat2_worse) %>%
  dplyr::mutate(resistance_level = "reslevel_1")
cat("  reslevel_1 (ros1 controls pat1, but not pat2):", nrow(df_reslevel_1), "parameter sets\n")

# reslevel_2: ros1_ controls pat2_ but NOT pat3_
df_step2 = df_step1 %>% dplyr::filter(!(param_set_id %in% df_reslevel_1$param_set_id))
df_reslevel_2 = df_step2 %>%
  dplyr::mutate(
    ros1_pat3_worse =       
      (mean_ros1_pat3_treg0_p>mean_ros1_pat2_treg0_p &
         mean_ros1_pat3_treg0_e<mean_ros1_pat2_treg0_e)
    
  ) %>%
  dplyr::filter(ros1_pat3_worse) %>%
  dplyr::mutate(resistance_level = "reslevel_2")
cat("  reslevel_2 (ros1 controls pat2, but not pat3):", nrow(df_reslevel_2), "parameter sets\n")

# ============= STEP 3: Identify evolutionary advantage of Tregs ====================================================
# 3.1: Within reslevel_1, tregs (treg1) should enable ros2_ to control pat2_ and/or pat3_
df_reslevel_1_with_treg_advantage = df_reslevel_1 %>%
  dplyr::mutate(
    ros2_treg1_controls_pat2 = (mean_ros2_pat2_treg0_p>mean_ros2_pat2_treg1_p & mean_ros2_pat2_treg0_e<mean_ros2_pat2_treg1_e),
    ros2_treg1_controls_pat3 = (mean_ros2_pat3_treg0_p>mean_ros2_pat3_treg1_p & mean_ros2_pat3_treg0_e<mean_ros2_pat3_treg1_e),
    # Tregs provide advantage if ros2+treg1 controls either pat2 or pat3 (or both)
    treg_advantage = ros2_treg1_controls_pat2 | ros2_treg1_controls_pat3
  ) %>%
  dplyr::filter(treg_advantage)

cat("  reslevel_1 with Treg advantage:", nrow(df_reslevel_1_with_treg_advantage), "parameter sets\n")

# 3.2: Within reslevel_2, tregs (treg1) should enable ros2_ to control pat3_
df_reslevel_2_with_treg_advantage = df_reslevel_2 %>%
  dplyr::mutate(
    ros2_treg1_controls_pat3 = (mean_ros2_pat3_treg0_p>mean_ros2_pat3_treg1_p & mean_ros2_pat3_treg0_e<mean_ros2_pat3_treg1_e),
    treg_advantage = ros2_treg1_controls_pat3
  ) %>%
  dplyr::filter(treg_advantage)

cat("  reslevel_2 with Treg advantage:", nrow(df_reslevel_2_with_treg_advantage), "parameter sets\n")

# ============= COMBINE AND SAVE RESULTS ====================================================
# Combine both resistance levels that show evolutionary selection for Tregs
df_evo_selected = bind_rows(
  df_reslevel_1_with_treg_advantage,
  df_reslevel_2_with_treg_advantage
)

cat("\nTotal evolutionarily selected parameter sets:", nrow(df_evo_selected), "\n")

# Save results
saveRDS(df_evo_selected, 'evo_selected.rds')

# Also save separate files for each resistance level for further analysis
saveRDS(df_reslevel_1_with_treg_advantage$param_set_id, 'evo_selected_reslevel_1.rds')
saveRDS(df_reslevel_2_with_treg_advantage$param_set_id, 'evo_selected_reslevel_2.rds')

# Print summary statistics
cat("\n=== SUMMARY ===\n")
cat("Resistance Level 1:", nrow(df_reslevel_1_with_treg_advantage), "sets\n")
cat("Resistance Level 2:", nrow(df_reslevel_2_with_treg_advantage), "sets\n")
cat("Total selected:", nrow(df_evo_selected), "sets\n")

df_reslevel_1_with_treg_advantage_long = df_reslevel_1_with_treg_advantage %>%
  dplyr::select(-resistance_level) %>%
  pivot_longer(
    cols = -c(param_set_id, injury_type),  # keep param_set_id as identifier
    names_to = "parameter",
    values_to = "value"
  )

df_params_pick = df_params %>% dplyr::filter(param_set_id %in% df_reslevel_1_with_treg_advantage$param_set_id)
df_params_pick_long = df_params_pick %>%
  pivot_longer(
    cols = -param_set_id,  # keep param_set_id as identifier
    names_to = "parameter",
    values_to = "value"
  )
