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
data_suffix    = '_ros_vs_ctrl' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
if(data_suffix == '_ros_vs_ctrl'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else{
  keep_param_id = reps_df %>% dplyr::filter(Freq==100) %>% dplyr::pull(Var1) # 100 = 10 reps per scenario, 5 scenarios x 2 times recording for epithelial and pathogen scores 
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
  mean_ctrl_pathogen_e,
  mean_tregs_off_pathogen_e,
  mean_ctrl_pathogen_p,
  mean_tregs_off_pathogen_p
))
df_comparisons_keep = df_comparisons

# ============= FILTER BASED ON CONTROL ====================================================
df_comparisons_ctrl_test = df_comparisons_keep %>% dplyr::filter(injury_type=='pathogenic') %>%
  dplyr::mutate(diff_ctrl_vs_tregs_off_p = mean_ctrl_pathogen_p-mean_tregs_off_pathogen_p,
                diff_ctrl_vs_tregs_off_e = mean_ctrl_pathogen_e-mean_tregs_off_pathogen_e)
df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test[c('param_set_id',
                                                             'diff_ctrl_vs_tregs_off_p',
                                                             'diff_ctrl_vs_tregs_off_e',
                                                             'mean_ctrl_pathogen_p',
                                                             'mean_ctrl_pathogen_e',
                                                             'mean_tregs_off_pathogen_p',
                                                             'mean_tregs_off_pathogen_e')]

df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple, 
                                        distinct(df_results[c('param_set_id',
                                                              'd_ctrl_vs_tregs_off_pathogen_p',
                                                              'd_ctrl_vs_tregs_off_pathogen_e')]), by='param_set_id')

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_better_p = ifelse(abs(d_ctrl_vs_tregs_off_pathogen_p)>=jsd_th 
                                    & diff_ctrl_vs_tregs_off_p>tol_in_p, 1, 
                                    ifelse(abs(d_ctrl_vs_tregs_off_pathogen_p)>=jsd_th 
                                           & diff_ctrl_vs_tregs_off_p < -1*tol_in_p,-1,0)))

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% 
  dplyr::mutate(ros_better_e = ifelse(abs(d_ctrl_vs_tregs_off_pathogen_e)>=jsd_th 
                                    & diff_ctrl_vs_tregs_off_e>tol_in_e, -1, 
                                    ifelse(abs(d_ctrl_vs_tregs_off_pathogen_e)>=jsd_th 
                                           & diff_ctrl_vs_tregs_off_e < -1*tol_in_e,1,0)))


df_comparisons_ctrl_test_simple_selected = df_comparisons_ctrl_test_simple %>% dplyr::filter(ros_better_e+ros_better_p==2)
df_comparisons_ctrl_test_simple_selected = df_comparisons_ctrl_test_simple_selected %>% dplyr::filter(mean_tregs_off_pathogen_e>149)
saveRDS(df_comparisons_ctrl_test_simple_selected,'evo_selected.rds')
dim(df_comparisons_ctrl_test_simple_selected)[1]
