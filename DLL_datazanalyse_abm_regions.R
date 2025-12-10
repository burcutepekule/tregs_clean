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
df_results_keep = readRDS('./data_cpp_read_abm.rds')
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
keep_param_id = reps_df %>% dplyr::filter(Freq==2000) %>% dplyr::pull(Var1) # 200 = 100 reps per scenario, 10 scenarios x 2 times recording for epithelial and pathogen scores 
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

source('./MISC/PLOT_REGIONS.R')
