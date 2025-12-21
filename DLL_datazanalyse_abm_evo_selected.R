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

df_comparisons = distinct(df_results %>% dplyr::select(
  param_set_id, injury_type,
  # Mean scores
  mean_ros0_pat1_treg0_e, # control
  mean_ros1_pat1_treg0_e, # ros on, pat level 1
  mean_ros1_pat2_treg0_e, # ros on, pat level 2
  mean_ros1_pat3_treg0_e, # ros on, pat level 3
  mean_ros2_pat1_treg0_e, # ros harder, pat level 1
  mean_ros2_pat2_treg0_e, # ros harder, pat level 2
  mean_ros2_pat3_treg0_e, # ros harder, pat level 3
  mean_ros1_pat1_treg1_e, # ros on, pat level 1, tregs on
  mean_ros1_pat2_treg1_e, # ros on, pat level 2, tregs on
  mean_ros1_pat3_treg1_e, # ros on, pat level 3, tregs on
  mean_ros2_pat1_treg1_e, # ros harder, pat level 1, tregs on
  mean_ros2_pat2_treg1_e, # ros harder, pat level 2, tregs on
  mean_ros2_pat3_treg1_e, # ros harder, pat level 3, tregs on
  mean_ros0_pat1_treg0_p, # control
  mean_ros1_pat1_treg0_p, # ros on, pat level 1
  mean_ros1_pat2_treg0_p, # ros on, pat level 2
  mean_ros1_pat3_treg0_p, # ros on, pat level 3
  mean_ros2_pat1_treg0_p, # ros harder, pat level 1
  mean_ros2_pat2_treg0_p, # ros harder, pat level 2
  mean_ros2_pat3_treg0_p, # ros harder, pat level 3
  mean_ros1_pat1_treg1_p, # ros on, pat level 1, tregs on
  mean_ros1_pat2_treg1_p, # ros on, pat level 2, tregs on
  mean_ros1_pat3_treg1_p, # ros on, pat level 3, tregs on
  mean_ros2_pat1_treg1_p, # ros harder, pat level 1, tregs on
  mean_ros2_pat2_treg1_p, # ros harder, pat level 2, tregs on
  mean_ros2_pat3_treg1_p  # ros harder, pat level 3, tregs on
))

df_comparisons_keep = df_comparisons

# ============= FILTER BASED ON CONTROL ====================================================
df_comparisons_p = df_comparisons_keep %>% dplyr::filter(injury_type=='pathogenic') %>%
  dplyr::mutate(diff_ros0_vs_ros1_pat1_tregs_off_p = mean_ros0_pat1_treg0_p-mean_ros1_pat1_treg0_p, # if this is higher, more pathogens when no ROS, thus ROS is BETTER
                diff_ros1_pat2_vs_pat1_tregs_off_p = mean_ros1_pat2_treg0_p-mean_ros1_pat1_treg0_p, # if this is higher, ROS_1 level is not enough for pat 2 level
                diff_ros2_pat2_vs_pat1_tregs_off_p = mean_ros2_pat2_treg0_p-mean_ros1_pat1_treg0_p,
                diff_ros1_pat3_vs_pat1_tregs_off_p = mean_ros1_pat3_treg0_p-mean_ros1_pat1_treg0_p, 
                diff_ros2_pat3_vs_pat1_tregs_off_p = mean_ros2_pat3_treg0_p-mean_ros1_pat1_treg0_p,
                diff_ros1_pat3_vs_pat2_tregs_off_p = mean_ros1_pat3_treg0_p-mean_ros1_pat2_treg0_p, 
                diff_ros2_pat3_vs_pat2_tregs_off_p = mean_ros2_pat3_treg0_p-mean_ros1_pat2_treg0_p
                )


df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test[c('param_set_id',
                                                             'diff_ros0_vs_ros1_pat1_tregs_off_p',
                                                             'diff_ros0_vs_ros1_pat1_tregs_off_e',
                                                             'mean_ros0_pat1_treg0_p',
                                                             'mean_ros0_pat1_treg0_e',
                                                             'mean_ros1_pat1_treg0_p',
                                                             'mean_ros1_pat1_treg0_e')]

df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple, 
                                        distinct(df_results[c('param_set_id',
                                                              'd_ros0_pat1_treg0_vs_ros1_pat1_treg0_p',
                                                              'd_ros0_pat1_treg0_vs_ros1_pat1_treg0_e')]), by='param_set_id')


### what is the idea here? ROS should be needed for the LOWEST pathogen level, that's for sure. But maybe that amount of ROS 
### is not enough for the second or third level?

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_better_p = ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_p)>=jsd_th 
                                    & diff_ros0_vs_ros1_pat1_tregs_off_p>tol_in_p, 1, 
                                    ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_p)>=jsd_th 
                                           & diff_ros0_vs_ros1_pat1_tregs_off_p < -1*tol_in_p,-1,0)))

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_better_e = ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_e)>=jsd_th 
                                    & diff_ros0_vs_ros1_pat1_tregs_off_e>tol_in_e, -1, 
                                    ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_e)>=jsd_th 
                                           & diff_ros0_vs_ros1_pat1_tregs_off_e < -1*tol_in_e,1,0)))

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_not_enough_pat2_p = ifelse(abs(d_ros1_pat1_treg0_vs_ros1_pat2_treg0_p)>=jsd_th 
                                      & diff_ros0_vs_ros1_pat1_tregs_off_p>tol_in_p, 1, 
                                      ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_p)>=jsd_th 
                                             & diff_ros0_vs_ros1_pat1_tregs_off_p < -1*tol_in_p,-1,0)))

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_not_enough_pat2_e = ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_e)>=jsd_th 
                                      & diff_ros0_vs_ros1_pat1_tregs_off_e>tol_in_e, -1, 
                                      ifelse(abs(d_ros0_pat1_treg0_vs_ros1_pat1_treg0_e)>=jsd_th 
                                             & diff_ros0_vs_ros1_pat1_tregs_off_e < -1*tol_in_e,1,0)))

df_comparisons_ctrl_test_simple_selected = df_comparisons_ctrl_test_simple %>% dplyr::filter(ros_better_e+ros_better_p==2)
df_comparisons_ctrl_test_simple_selected = df_comparisons_ctrl_test_simple_selected %>% dplyr::filter(mean_ros1_pat1_treg0_e>149)
saveRDS(df_comparisons_ctrl_test_simple_selected,'evo_selected.rds')
dim(df_comparisons_ctrl_test_simple_selected)[1]
