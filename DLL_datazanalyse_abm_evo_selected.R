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

### ======== Any treg help for other patros levels?
df_step1_treg_jsd = df_step1[c('param_set_id',
                           'd_ros1_pat1_treg0_vs_ros1_pat1_treg1_e',
                           'd_ros2_pat1_treg0_vs_ros2_pat1_treg1_e',
                           'd_ros1_pat2_treg0_vs_ros1_pat2_treg1_e',
                           'd_ros2_pat2_treg0_vs_ros2_pat2_treg1_e',
                           'd_ros1_pat3_treg0_vs_ros1_pat3_treg1_e',
                           'd_ros2_pat3_treg0_vs_ros2_pat3_treg1_e',
                           'd_ros1_pat1_treg0_vs_ros1_pat1_treg1_p',
                           'd_ros2_pat1_treg0_vs_ros2_pat1_treg1_p',
                           'd_ros1_pat2_treg0_vs_ros1_pat2_treg1_p',
                           'd_ros2_pat2_treg0_vs_ros2_pat2_treg1_p',
                           'd_ros1_pat3_treg0_vs_ros1_pat3_treg1_p',
                           'd_ros2_pat3_treg0_vs_ros2_pat3_treg1_p')]

df_step1_treg_mean = df_step1[c('param_set_id',
                               'mean_ros1_pat1_treg0_e','mean_ros1_pat1_treg1_e',
                               'mean_ros2_pat1_treg0_e','mean_ros2_pat1_treg1_e',
                               'mean_ros1_pat2_treg0_e','mean_ros1_pat2_treg1_e',
                               'mean_ros2_pat2_treg0_e','mean_ros2_pat2_treg1_e',
                               'mean_ros1_pat3_treg0_e','mean_ros1_pat3_treg1_e',
                               'mean_ros2_pat3_treg0_e','mean_ros2_pat3_treg1_e',
                               'mean_ros1_pat1_treg0_p','mean_ros1_pat1_treg1_p',
                               'mean_ros2_pat1_treg0_p','mean_ros2_pat1_treg1_p',
                               'mean_ros1_pat2_treg0_p','mean_ros1_pat2_treg1_p',
                               'mean_ros2_pat2_treg0_p','mean_ros2_pat2_treg1_p',
                               'mean_ros1_pat3_treg0_p','mean_ros1_pat3_treg1_p',
                               'mean_ros2_pat3_treg0_p','mean_ros2_pat3_treg1_p')]

df_step1_treg_mean = df_step1_treg_mean %>% dplyr::rowwise() %>% dplyr::mutate(
  diff_ros1_pat1_e = mean_ros1_pat1_treg1_e-mean_ros1_pat1_treg0_e,
  diff_ros2_pat1_e = mean_ros2_pat1_treg1_e-mean_ros2_pat1_treg0_e,
  diff_ros1_pat2_e = mean_ros1_pat2_treg1_e-mean_ros1_pat2_treg0_e,
  diff_ros2_pat2_e = mean_ros2_pat2_treg1_e-mean_ros2_pat2_treg0_e,
  diff_ros1_pat3_e = mean_ros1_pat3_treg1_e-mean_ros1_pat3_treg0_e,
  diff_ros2_pat3_e = mean_ros2_pat3_treg1_e-mean_ros2_pat3_treg0_e,
  diff_ros1_pat1_p = mean_ros1_pat1_treg1_p-mean_ros1_pat1_treg0_p,
  diff_ros2_pat1_p = mean_ros2_pat1_treg1_p-mean_ros2_pat1_treg0_p,
  diff_ros1_pat2_p = mean_ros1_pat2_treg1_p-mean_ros1_pat2_treg0_p,
  diff_ros2_pat2_p = mean_ros2_pat2_treg1_p-mean_ros2_pat2_treg0_p,
  diff_ros1_pat3_p = mean_ros1_pat3_treg1_p-mean_ros1_pat3_treg0_p,
  diff_ros2_pat3_p = mean_ros2_pat3_treg1_p-mean_ros2_pat3_treg0_p
)

df_step1_treg_merged = merge(df_step1_treg_jsd, df_step1_treg_mean[c('param_set_id',
                                              'diff_ros1_pat1_e',
                                              'diff_ros2_pat1_e',
                                              'diff_ros1_pat2_e',
                                              'diff_ros2_pat2_e',
                                              'diff_ros1_pat3_e',
                                              'diff_ros2_pat3_e',
                                              'diff_ros1_pat1_p',
                                              'diff_ros2_pat1_p',
                                              'diff_ros1_pat2_p',
                                              'diff_ros2_pat2_p',
                                              'diff_ros1_pat3_p',
                                              'diff_ros2_pat3_p')], by='param_set_id')


df_step1_treg_merged = df_step1_treg_merged %>% dplyr::mutate(
  ros1_pat1_e = ifelse(diff_ros1_pat1_e>tol_in_e & d_ros1_pat1_treg0_vs_ros1_pat1_treg1_e>jsd_th, T, F),
  ros2_pat1_e = ifelse(diff_ros2_pat1_e>tol_in_e & d_ros2_pat1_treg0_vs_ros2_pat1_treg1_e>jsd_th, T, F),
  ros1_pat2_e = ifelse(diff_ros1_pat2_e>tol_in_e & d_ros1_pat2_treg0_vs_ros1_pat2_treg1_e>jsd_th, T, F),
  ros2_pat2_e = ifelse(diff_ros2_pat2_e>tol_in_e & d_ros2_pat2_treg0_vs_ros2_pat2_treg1_e>jsd_th, T, F),
  ros1_pat3_e = ifelse(diff_ros1_pat3_e>tol_in_e & d_ros1_pat3_treg0_vs_ros1_pat3_treg1_e>jsd_th, T, F),
  ros2_pat3_e = ifelse(diff_ros2_pat3_e>tol_in_e & d_ros2_pat3_treg0_vs_ros2_pat3_treg1_e>jsd_th, T, F),
  ros1_pat1_p = ifelse(diff_ros1_pat1_p<tol_in_p & d_ros1_pat1_treg0_vs_ros1_pat1_treg1_p>jsd_th, T, F),
  ros2_pat1_p = ifelse(diff_ros2_pat1_p<tol_in_p & d_ros2_pat1_treg0_vs_ros2_pat1_treg1_p>jsd_th, T, F),
  ros1_pat2_p = ifelse(diff_ros1_pat2_p<tol_in_p & d_ros1_pat2_treg0_vs_ros1_pat2_treg1_p>jsd_th, T, F),
  ros2_pat2_p = ifelse(diff_ros2_pat2_p<tol_in_p & d_ros2_pat2_treg0_vs_ros2_pat2_treg1_p>jsd_th, T, F),
  ros1_pat3_p = ifelse(diff_ros1_pat3_p<tol_in_p & d_ros1_pat3_treg0_vs_ros1_pat3_treg1_p>jsd_th, T, F),
  ros2_pat3_p = ifelse(diff_ros2_pat3_p<tol_in_p & d_ros2_pat3_treg0_vs_ros2_pat3_treg1_p>jsd_th, T, F)
)

df_step1_treg_merged_slim = df_step1_treg_merged[c('param_set_id',
                                                   'ros1_pat1_e','ros2_pat1_e',
                                                   'ros1_pat2_e','ros2_pat2_e',
                                                   'ros1_pat3_e','ros2_pat3_e',
                                                   'ros1_pat1_p','ros2_pat1_p',
                                                   'ros1_pat2_p','ros2_pat2_p',
                                                   'ros1_pat3_p','ros2_pat3_p')]
# at least one TRUE in them?
df_with_true = df_step1_treg_merged_slim[apply(df_step1_treg_merged_slim[, -1], 1, any), ]

saveRDS(df_with_true$param_set_id, 'evo_selected_js_based.rds')
dim(df_with_true)

