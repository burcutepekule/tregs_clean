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
tol_in_p       = -1*25*25*0.05
# tol_in_e       = 0
# tol_in_p       = 0
M1_M2_diff     = 0
filter_control = 0
labels_on      = 1
score_type     = 'pathogen'
data_suffix    = '_ros_vs_ctrl_nopatlevel' #
data_suffix    = '_ros_vs_ctrl_patros' #
data_suffix    = '_patros' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
unique(df_results_keep$param_set_id)
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df = as.data.frame(table(df_results_keep$param_set_id))
reps_df$Var1 = as.numeric(as.character(reps_df$Var1))
if(data_suffix == '_ros_vs_ctrl_nopatlevel'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else if(data_suffix == '_ros_vs_ctrl_patros'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==260) %>% dplyr::pull(Var1) # 260 = 10 reps per scenario, 13 scenarios x 2 times recording for epithelial and pathogen scores 
}else if(data_suffix == '_patros'){
  # keep_param_id = reps_df %>% dplyr::filter(Freq==1100) %>% dplyr::pull(Var1) # 1100 = 5 reps per scenario, 110 scenarios x 2 times recording for epithelial and pathogen scores
  keep_param_id = reps_df %>% dplyr::filter(Freq > 440) %>% dplyr::pull(Var1) # 1100 = 5 reps per scenario, 110 scenarios x 2 times recording for epithelial and pathogen scores
}

df_results = df_results_keep %>% filter(param_set_id %in% keep_param_id)
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
# ss_start_threshold = 4500 # used to be 4500, just for simulation purposes to save time
# ss_start_threshold = 1800 # used to be 4500, just for simulation purposes to save time
# param_id_all_below = df_results %>%
#   dplyr::group_by(param_set_id) %>%
#   dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
#   dplyr::filter(all_below) %>%
#   dplyr::pull(param_set_id)
# df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
# num_params = length(unique(df_results$param_set_id))
# print(num_params)
# saveRDS(unique(df_results$param_set_id), 'evo_selected_level_0.rds')

df_comparisons = df_results %>% dplyr::select(
  param_set_id, injury_type, 
  starts_with("d_ros"), starts_with("mean_ros"))

df_comparisons_keep = distinct(df_comparisons)

df_comparisons_keep_pick = df_comparisons_keep %>% dplyr::filter(param_set_id==4201) %>% dplyr::select(-injury_type)

df_comparisons_keep_pick = pivot_longer(
  df_comparisons_keep_pick, 
  cols = -c(param_set_id),  # or specify which cols TO pivot
  names_to = 'variable',
  values_to = 'value'
)
# ============= HELPER FUNCTION: Define "infection under control" ====================================================
# Infection is under control when:
# 1. Epithelial health is approximately 150 (threshold: > 149)
# 2. Pathogen load is approximately 0 (threshold: < 1)
is_under_control = function(epithelial_health, pathogen_load) {
  return(epithelial_health > 149 & pathogen_load < 2)
}

# ============= STEP 1: Filter for evolutionary need for ROS ====================================================
# Rule 1: Without ROS (ros_0, treg0_), even the lowest pathogen level (pat1_) should:
#   a) Cause significant injury to epithelium (ros0_pat1 NOT under control)
#   b) Have higher pathogen burden than ros1_pat1 case (ros1_pat1 should be better)

# Start with the base filtering
df_steps = df_comparisons_keep %>%
  dplyr::filter(injury_type == 'pathogenic')

ros_vals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
pat_vals = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Dynamically create the controlled status columns for all ros/pat combinations
for (ros in ros_vals) {
  for (pat in pat_vals) {
    col_name = paste0("ros", ros, "_pat", pat, "_controlled")
    mean_e_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_e")
    mean_p_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_p")
    
    df_steps = df_steps %>%
      dplyr::mutate(!!col_name := is_under_control(.data[[mean_e_col]], .data[[mean_p_col]]))
  }
}


# For each ros level, check if the FALSE -> TRUE -> FALSE pattern holds
# This adds a new column for each ros level indicating if the pattern is present

for (ros in ros_vals) {
  # Create column names for this ros level across all pat levels
  controlled_cols = paste0("ros", ros, "_pat", pat_vals, "_controlled")
  
  # Create a new column to store whether the pattern holds for this ros level
  pattern_col_name = paste0("ros", ros, "_pattern_FTF")
  
  df_steps[[pattern_col_name]] = apply(df_steps[, controlled_cols], 1, function(row) {
    # Find indices where value is TRUE
    true_indices = which(row)
    
    if (length(true_indices) == 0) {
      # No TRUE values, pattern doesn't hold
      return(FALSE)
    }
    
    first_true = min(true_indices)
    last_true = max(true_indices)
    
    # Pattern holds if:
    # 1. There's at least one FALSE before the first TRUE (first_true > 1)
    # 2. There's at least one FALSE after the last TRUE (last_true < length)
    # 3. There's at least one TRUE (already checked above)
    pattern_holds = (first_true > 1) && (last_true < length(row))
    
    return(pattern_holds)
  })
}

# Now you can filter to find parameter sets that show this pattern
# For example, to see which param sets have the pattern for ros level 0:
df_with_pattern = df_steps %>% dplyr::filter(ros3_pattern_FTF == TRUE)

# Or to see how many param sets have the pattern for each ros level:
pattern_summary = data.frame(
  ros_level = ros_vals,
  n_with_pattern = sapply(ros_vals, function(ros) {
    sum(df_steps[[paste0("ros", ros, "_pattern_FTF")]], na.rm = TRUE)
  })
)
print(pattern_summary)









# 
# length(df_steps_not_selected$param_set_id)
# 
# df_steps_not_selected_0 = df_steps_not_selected %>% filter(ros0_pat1_controlled==T)
# 
# df_step_pat1 = df_steps %>% dplyr::filter(!ros0_pat1_controlled)
# df_step_pat2 = df_steps %>% dplyr::filter(ros0_pat1_controlled & !ros0_pat2_controlled)
# df_step_pat3 = df_steps %>% dplyr::filter(ros0_pat1_controlled & ros0_pat2_controlled & !ros0_pat3_controlled)
# 
# length(df_comparisons_keep$param_set_id)
# length(df_step_pat1$param_set_id) # without ROS, no level can be controlled 
# # - so minimum pathogen load to take seriously for these parameter vectors is pat level 1
# # can tregs help control pathogen level 1 to begin with? and then 2 and/or 3?
# length(df_step_pat2$param_set_id) # without ROS, only pat level 1 can be controlled, but not pat levels 2 and 3 
# # - so minimum pathogen load to take seriously for these parameter vectors is pat level 2
# # can tregs help control pathogen level 2 to begin with? and then 3?
# length(df_step_pat3$param_set_id) # without ROS, only pat levels 1 & 2 can be controlled, but not pat level 3
# # - so minimum pathogen load to take seriously for these parameter vectors is pat level 3
# # can tregs help control pathogen level 3?
# 
# evo_selected_level_1_ids = sort(df_step_pat1$param_set_id)
# evo_selected_level_2_ids = sort(df_step_pat2$param_set_id)
# evo_selected_level_3_ids = sort(df_step_pat3$param_set_id)
# 
# saveRDS(evo_selected_level_1_ids, 'evo_selected_level_1.rds')
# saveRDS(evo_selected_level_2_ids, 'evo_selected_level_2.rds')
# saveRDS(evo_selected_level_3_ids, 'evo_selected_level_3.rds')
# 
# ###### HERE - NOW I HAVE TO THINK ABOUT DIFFERENT LEVELS OF ROS WORKING OR NOT GIVEN DIFFERENT LEVELS OF PATHOGENS!
# 
# 
# ## OK, find cases where
# 
# ## Option 1 
# ## ROS 0 (control) cannot resolve the infection for pat level 1. 
# ## ROS 1 can resolve the infection for pat level 1. (Because at one point they needto be selected by evolution to continue before Tregs evolve.)
# ## ROS 1 cannot resolve the infection for pat level 2. That's where the organism fails.
# ## 1.a ROS 2 cannot resolve the infection for pat level 2 on its own but can resolve the infection with the help of Tregs.
# ## 1.b ROS 3 cannot resolve the infection for pat level 2 on its own but can resolve the infection with the help of Tregs.
# 
# ## Option 2 
# ## ROS 0 (control) cannot resolve the infection for pat level 2. 
# ## ROS 1 can resolve the infection for pat level 1. (Because at one point they needto be selected by evolution to continue before Tregs evolve.)
# ## ROS 1 cannot resolve the infection for pat level 2. That's where the organism fails.
# ## 1.a ROS 2 cannot resolve the infection on its own but can resolve the infection with the help of Tregs.
# ## 1.b ROS 3 cannot resolve the infection on its own but can resolve the infection with the help of Tregs.
# 
# 
# 
# 
# 
# # 
# # ### ======== Any treg help for other patros levels?
# # df_step1_treg_jsd = df_step1[c('param_set_id',
# #                            'd_ros1_pat1_treg0_vs_ros1_pat1_treg1_e',
# #                            'd_ros2_pat1_treg0_vs_ros2_pat1_treg1_e',
# #                            'd_ros1_pat2_treg0_vs_ros1_pat2_treg1_e',
# #                            'd_ros2_pat2_treg0_vs_ros2_pat2_treg1_e',
# #                            'd_ros1_pat3_treg0_vs_ros1_pat3_treg1_e',
# #                            'd_ros2_pat3_treg0_vs_ros2_pat3_treg1_e',
# #                            'd_ros1_pat1_treg0_vs_ros1_pat1_treg1_p',
# #                            'd_ros2_pat1_treg0_vs_ros2_pat1_treg1_p',
# #                            'd_ros1_pat2_treg0_vs_ros1_pat2_treg1_p',
# #                            'd_ros2_pat2_treg0_vs_ros2_pat2_treg1_p',
# #                            'd_ros1_pat3_treg0_vs_ros1_pat3_treg1_p',
# #                            'd_ros2_pat3_treg0_vs_ros2_pat3_treg1_p')]
# # 
# # df_step1_treg_mean = df_step1[c('param_set_id',
# #                                'mean_ros1_pat1_treg0_e','mean_ros1_pat1_treg1_e',
# #                                'mean_ros2_pat1_treg0_e','mean_ros2_pat1_treg1_e',
# #                                'mean_ros1_pat2_treg0_e','mean_ros1_pat2_treg1_e',
# #                                'mean_ros2_pat2_treg0_e','mean_ros2_pat2_treg1_e',
# #                                'mean_ros1_pat3_treg0_e','mean_ros1_pat3_treg1_e',
# #                                'mean_ros2_pat3_treg0_e','mean_ros2_pat3_treg1_e',
# #                                'mean_ros1_pat1_treg0_p','mean_ros1_pat1_treg1_p',
# #                                'mean_ros2_pat1_treg0_p','mean_ros2_pat1_treg1_p',
# #                                'mean_ros1_pat2_treg0_p','mean_ros1_pat2_treg1_p',
# #                                'mean_ros2_pat2_treg0_p','mean_ros2_pat2_treg1_p',
# #                                'mean_ros1_pat3_treg0_p','mean_ros1_pat3_treg1_p',
# #                                'mean_ros2_pat3_treg0_p','mean_ros2_pat3_treg1_p')]
# # 
# # df_step1_treg_mean = df_step1_treg_mean %>% dplyr::rowwise() %>% dplyr::mutate(
# #   diff_ros1_pat1_e = mean_ros1_pat1_treg1_e-mean_ros1_pat1_treg0_e,
# #   diff_ros2_pat1_e = mean_ros2_pat1_treg1_e-mean_ros2_pat1_treg0_e,
# #   diff_ros1_pat2_e = mean_ros1_pat2_treg1_e-mean_ros1_pat2_treg0_e,
# #   diff_ros2_pat2_e = mean_ros2_pat2_treg1_e-mean_ros2_pat2_treg0_e,
# #   diff_ros1_pat3_e = mean_ros1_pat3_treg1_e-mean_ros1_pat3_treg0_e,
# #   diff_ros2_pat3_e = mean_ros2_pat3_treg1_e-mean_ros2_pat3_treg0_e,
# #   diff_ros1_pat1_p = mean_ros1_pat1_treg1_p-mean_ros1_pat1_treg0_p,
# #   diff_ros2_pat1_p = mean_ros2_pat1_treg1_p-mean_ros2_pat1_treg0_p,
# #   diff_ros1_pat2_p = mean_ros1_pat2_treg1_p-mean_ros1_pat2_treg0_p,
# #   diff_ros2_pat2_p = mean_ros2_pat2_treg1_p-mean_ros2_pat2_treg0_p,
# #   diff_ros1_pat3_p = mean_ros1_pat3_treg1_p-mean_ros1_pat3_treg0_p,
# #   diff_ros2_pat3_p = mean_ros2_pat3_treg1_p-mean_ros2_pat3_treg0_p
# # )
# # 
# # df_step1_treg_merged = merge(df_step1_treg_jsd, df_step1_treg_mean[c('param_set_id',
# #                                               'diff_ros1_pat1_e',
# #                                               'diff_ros2_pat1_e',
# #                                               'diff_ros1_pat2_e',
# #                                               'diff_ros2_pat2_e',
# #                                               'diff_ros1_pat3_e',
# #                                               'diff_ros2_pat3_e',
# #                                               'diff_ros1_pat1_p',
# #                                               'diff_ros2_pat1_p',
# #                                               'diff_ros1_pat2_p',
# #                                               'diff_ros2_pat2_p',
# #                                               'diff_ros1_pat3_p',
# #                                               'diff_ros2_pat3_p')], by='param_set_id')
# # 
# # 
# # df_step1_treg_merged = df_step1_treg_merged %>% dplyr::mutate(
# #   ros1_pat1_e = ifelse(diff_ros1_pat1_e>tol_in_e & d_ros1_pat1_treg0_vs_ros1_pat1_treg1_e>jsd_th, T, F),
# #   ros2_pat1_e = ifelse(diff_ros2_pat1_e>tol_in_e & d_ros2_pat1_treg0_vs_ros2_pat1_treg1_e>jsd_th, T, F),
# #   ros1_pat2_e = ifelse(diff_ros1_pat2_e>tol_in_e & d_ros1_pat2_treg0_vs_ros1_pat2_treg1_e>jsd_th, T, F),
# #   ros2_pat2_e = ifelse(diff_ros2_pat2_e>tol_in_e & d_ros2_pat2_treg0_vs_ros2_pat2_treg1_e>jsd_th, T, F),
# #   ros1_pat3_e = ifelse(diff_ros1_pat3_e>tol_in_e & d_ros1_pat3_treg0_vs_ros1_pat3_treg1_e>jsd_th, T, F),
# #   ros2_pat3_e = ifelse(diff_ros2_pat3_e>tol_in_e & d_ros2_pat3_treg0_vs_ros2_pat3_treg1_e>jsd_th, T, F),
# #   ros1_pat1_p = ifelse(diff_ros1_pat1_p<tol_in_p & d_ros1_pat1_treg0_vs_ros1_pat1_treg1_p>jsd_th, T, F),
# #   ros2_pat1_p = ifelse(diff_ros2_pat1_p<tol_in_p & d_ros2_pat1_treg0_vs_ros2_pat1_treg1_p>jsd_th, T, F),
# #   ros1_pat2_p = ifelse(diff_ros1_pat2_p<tol_in_p & d_ros1_pat2_treg0_vs_ros1_pat2_treg1_p>jsd_th, T, F),
# #   ros2_pat2_p = ifelse(diff_ros2_pat2_p<tol_in_p & d_ros2_pat2_treg0_vs_ros2_pat2_treg1_p>jsd_th, T, F),
# #   ros1_pat3_p = ifelse(diff_ros1_pat3_p<tol_in_p & d_ros1_pat3_treg0_vs_ros1_pat3_treg1_p>jsd_th, T, F),
# #   ros2_pat3_p = ifelse(diff_ros2_pat3_p<tol_in_p & d_ros2_pat3_treg0_vs_ros2_pat3_treg1_p>jsd_th, T, F)
# # )
# # 
# # df_step1_treg_merged_slim = df_step1_treg_merged[c('param_set_id',
# #                                                    'ros1_pat1_e','ros2_pat1_e',
# #                                                    'ros1_pat2_e','ros2_pat2_e',
# #                                                    'ros1_pat3_e','ros2_pat3_e',
# #                                                    'ros1_pat1_p','ros2_pat1_p',
# #                                                    'ros1_pat2_p','ros2_pat2_p',
# #                                                    'ros1_pat3_p','ros2_pat3_p')]
# # # at least one TRUE in them?
# # df_with_true = df_step1_treg_merged_slim[apply(df_step1_treg_merged_slim[, -1], 1, any), ]
# # 
# # saveRDS(df_with_true$param_set_id, 'evo_selected_js_based.rds')
# # dim(df_with_true)
# # 
