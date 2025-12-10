### FINDS CASES WHERE EPITHELIAL SCORE IS BETTER/WORSE AND PATHOGENIC LOAD IS WORSE/BETTER
rm(list=ls())
jsd_th         = 0.3
tol_in         = 125*0.25
M1_M2_diff     = 1
filter_control = 1
labels_on      = 1
th_ctrl_pat    = -Inf
data_suffix    = '_10' # empty for 100 reps, _10 for 10 reps

analysis_pick  = 2

if(analysis_pick==1){
  # 1 =========================================
  condition_subt_from = 'ctrl'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'ctrl_vs_tregs_off'
}else if(analysis_pick==2){
  # 2 =========================================
  condition_subt_from = 'tregs_on'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'tregs_on_vs_off'
}else if(analysis_pick==3){
  # 3 =========================================
  condition_subt_from = 'tregs_on'
  condition_subt      = 'tregs_rnd'
  jensen_distance     = 'tregs_on_vs_rnd'
}else if(analysis_pick==4){
  # 4 =========================================
  condition_subt_from = 'macspec'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'macspec_vs_tregs_off'
}else if(analysis_pick==5){
  # 5 =========================================
  condition_subt_from = 'macspec'
  condition_subt      = 'tregs_on'
  jensen_distance     = 'macspec_vs_tregs_on'
}else if(analysis_pick==6){
  # 6 =========================================
  condition_subt_from = 'macspec'
  condition_subt      = 'tregs_rnd'
  jensen_distance     = 'macspec_vs_tregs_rnd'
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
if(data_suffix == '_10'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==200) %>% dplyr::pull(Var1) # 20 = 10 reps per scenario, 10 scenarios x 2 times recording for epithelial and pathogen scores 
}else{
  keep_param_id = reps_df %>% dplyr::filter(Freq==2000) %>% dplyr::pull(Var1) # 200 = 100 reps per scenario, 10 scenarios x 2 times recording for epithelial and pathogen scores 
}
df_results    = df_results_keep %>% filter(param_set_id %in% keep_param_id)

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id)) # >99%!
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold # TRUE, sanity check
unique(table(df_results$param_set_id)) #2000, sanity check
length(unique(df_results$param_set_id))

df_comparisons = distinct(df_results %>% dplyr::select(
  param_set_id, injury_type,
  # Mean scores
  mean_ctrl_sterile_e, mean_ctrl_pathogen_e,
  mean_tregs_off_sterile_e, mean_tregs_off_pathogen_e,
  mean_tregs_on_sterile_e, mean_tregs_on_pathogen_e,
  mean_tregs_rnd_sterile_e, mean_tregs_rnd_pathogen_e,
  mean_macspec_sterile_e, mean_macspec_pathogen_e,
  mean_ctrl_sterile_p, mean_ctrl_pathogen_p,
  mean_tregs_off_sterile_p, mean_tregs_off_pathogen_p,
  mean_tregs_on_sterile_p, mean_tregs_on_pathogen_p,
  mean_tregs_rnd_sterile_p, mean_tregs_rnd_pathogen_p,
  mean_macspec_sterile_p, mean_macspec_pathogen_p
))
df_comparisons_keep = df_comparisons

# ============= FILTER BASED ON CONTROL ====================================================
if(filter_control==1){ # This is to pick cases where ROS is an evolutionary favorable thing to have!
  df_comparisons_ctrl_test = df_comparisons_keep %>%
    dplyr::mutate(diff_ctrl_vs_tregs_off = mean_ctrl_pathogen_p-mean_tregs_off_pathogen_p)
  df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test[c('param_set_id','diff_ctrl_vs_tregs_off','mean_ctrl_pathogen_p','mean_tregs_off_pathogen_p')]
  df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple, distinct(df_results[c('param_set_id','d_ctrl_vs_tregs_off_pathogen_p')]), by='param_set_id')
  df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% dplyr::filter(abs(d_ctrl_vs_tregs_off_pathogen_p)>=jsd_th 
                                                                                      & diff_ctrl_vs_tregs_off>tol_in) # has to be positive, control should have more pathogens
  #& mean_ctrl_pathogen_p>th_ctrl_pat) # this is not necessary and also a bit arbitrary
  ids_matter     = unique(df_comparisons_ctrl_test_simple %>% dplyr::pull(param_set_id))
  df_comparisons = df_comparisons %>% dplyr::filter(param_set_id %in% ids_matter)
  # plot(df_comparisons_ctrl_test_simple$diff_ctrl_vs_tregs_off, df_comparisons_ctrl_test_simple$mean_ctrl_pathogen_p)
  # plot(df_comparisons_ctrl_test_simple$mean_tregs_off_pathogen_p, df_comparisons_ctrl_test_simple$mean_ctrl_pathogen_p)
  
}
# ============= FILTER BASED ON CONTROL ====================================================
df_comparisons_keep = df_comparisons
# ==========================================================================================
score_type = 'epithelial'
mean_subt_from_full_p = paste0('mean_',condition_subt_from,'_pathogen_',substr(score_type,1,1))
mean_subt_from_full_s = paste0('mean_',condition_subt_from,'_sterile_',substr(score_type,1,1))

mean_subt_full_p = paste0('mean_',condition_subt,'_pathogen_',substr(score_type,1,1))
mean_subt_full_s = paste0('mean_',condition_subt,'_sterile_',substr(score_type,1,1))

d_full_p = paste0('d_',jensen_distance,'_pathogen_',substr(score_type,1,1))
d_full_s = paste0('d_',jensen_distance,'_sterile_',substr(score_type,1,1))

df_comparisons = df_comparisons %>% 
  dplyr::mutate(diff_compare = ifelse(
    injury_type == 'pathogenic', 
    .data[[mean_subt_from_full_p]] - .data[[mean_subt_full_p]],
    .data[[mean_subt_from_full_s]] - .data[[mean_subt_full_s]]
  ))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_compare)
df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id',d_full_p, d_full_s)]), by='param_set_id')

# ====================== Find significant cases ======================================================
df_comparisons_plot   = df_comparisons
df_comparisons_plot_1 = df_comparisons_plot # just to keep them separate
# ====================== PATHOGENIC ==================================================================
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in, 1,
                                                                                                       ifelse(diff_compare < -1*tol_in,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                              diff_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'

# ====================== STERILE ==================================================================
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in, 1,
                                                                                                 ifelse(diff_compare < -1*tol_in,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_s]])>jsd_th,
                                                                                                        diff_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-all_of(d_full_p))
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)==d_full_s)]='cohens_d'

diff_better_sterile_e_score    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
diff_worse_sterile_e_score     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_better_pathogenic_e_score = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
diff_worse_pathogenic_e_score  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)

# ==========================================================================================
df_comparisons = df_comparisons_keep
score_type = 'pathogen'
# ==========================================================================================
mean_subt_from_full_p = paste0('mean_',condition_subt_from,'_pathogen_',substr(score_type,1,1))
mean_subt_from_full_s = paste0('mean_',condition_subt_from,'_sterile_',substr(score_type,1,1))

mean_subt_full_p = paste0('mean_',condition_subt,'_pathogen_',substr(score_type,1,1))
mean_subt_full_s = paste0('mean_',condition_subt,'_sterile_',substr(score_type,1,1))

d_full_p = paste0('d_',jensen_distance,'_pathogen_',substr(score_type,1,1))
d_full_s = paste0('d_',jensen_distance,'_sterile_',substr(score_type,1,1))

df_comparisons = df_comparisons %>% 
  dplyr::mutate(diff_compare = ifelse(
    injury_type == 'pathogenic', 
    .data[[mean_subt_from_full_p]] - .data[[mean_subt_full_p]],
    .data[[mean_subt_from_full_s]] - .data[[mean_subt_full_s]]
  ))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_compare)
df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id',d_full_p, d_full_s)]), by='param_set_id')

# ====================== Find significant cases ======================================================
df_comparisons_plot   = df_comparisons
df_comparisons_plot_2 = df_comparisons_plot # just to keep them separate
# ====================== PATHOGENIC ==================================================================
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in, 1,
                                                                                                       ifelse(diff_compare < -1*tol_in,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                              diff_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'

# ====================== STERILE ==================================================================
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in, 1,
                                                                                                 ifelse(diff_compare < -1*tol_in,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_s]])>jsd_th,
                                                                                                        diff_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-all_of(d_full_p))
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)==d_full_s)]='cohens_d'

# HERE -1 +1 SHOULD SWITCH BECAUSE HIGHER PATHOGEN MEANS WORSE FOR PATHOGEN CONTOL! IT'S NOT LIKE EPITHELIAL SCORE! 
diff_better_sterile_p_score    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_worse_sterile_p_score     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==+1) %>% dplyr::pull(param_set_id)
diff_better_pathogenic_p_score = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_worse_pathogenic_p_score  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==+1) %>% dplyr::pull(param_set_id)

# ============ CONFLICT ===========================================================================
# conflict can happen only for pathogenic cases since sterile cases do not have pathogens!
### BOTH ZERO!
intersect(diff_better_pathogenic_p_score, diff_worse_pathogenic_e_score) 
intersect(diff_worse_pathogenic_p_score, diff_better_pathogenic_e_score) 

### BETTER OR WORSE FOR BOTH CASES?
intersect(diff_better_pathogenic_p_score, diff_better_pathogenic_e_score) 
intersect(diff_worse_pathogenic_p_score, diff_worse_pathogenic_e_score) 

