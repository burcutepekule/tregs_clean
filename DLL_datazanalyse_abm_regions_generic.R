# Generic script for comparing two groups in ABM regions analysis
#
# Required parameters (must be set before sourcing this script):
# - comparison_name: Name of the comparison (e.g., "tregs_on_vs_tregs_off")
# - group1: Name of first group (e.g., "tregs_on")
# - group2: Name of second group (e.g., "tregs_off")
# - jsd_th: Jensen-Shannon distance threshold
# - tol_in: Tolerance for difference threshold (can be NA for auto-calculation)
# - filter_control: Whether to filter based on control (0 or 1)
#
# Optional parameters:
# - perform_randomness_check: Whether to check randomized tregs (default: FALSE)
# - perform_macspec_analysis: Whether to perform macspec-specific analysis (default: FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Check required parameters
required_params = c("comparison_name", "group1", "group2", "jsd_th", "filter_control")
missing_params = setdiff(required_params, ls())
if(length(missing_params) > 0) {
  stop(paste("Missing required parameters:", paste(missing_params, collapse=", ")))
}

# Set defaults for optional parameters
if(!exists("perform_randomness_check")) perform_randomness_check = FALSE
if(!exists("perform_macspec_analysis")) perform_macspec_analysis = FALSE
if(!exists("tol_in")) tol_in = NA

# ============= DATA LOADING AND FILTERING =============
df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS('./data_cpp_read_abm.rds')
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
keep_param_id = reps_df %>% dplyr::filter(Freq==200) %>% dplyr::pull(Var1)
df_results    = df_results_keep %>% filter(param_set_id %in% keep_param_id)
unique(table(df_results$param_set_id))
length(unique(df_results$param_set_id))

#----- filter based on ss_start
ss_start_threshold = 4500
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id))
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold
unique(table(df_results$param_set_id))
length(unique(df_results$param_set_id))

# ============= PREPARE COMPARISON DATA =============
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

# ============= FILTER BASED ON CONTROL =============
if(filter_control==1){
  df_comparisons_ctrl_test = df_comparisons_keep %>%
    dplyr::mutate(diff_ctrl_test = mean_ctrl_pathogen_p-mean_tregs_off_pathogen_p)

  df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test[c('param_set_id','diff_ctrl_test')]
  df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple,
                                          distinct(df_results[c('param_set_id','d_ctrl_vs_test_pathogen_p')]),
                                          by='param_set_id')
  ids_matter_df = df_comparisons_ctrl_test_simple %>%
    dplyr::filter(abs(d_ctrl_vs_test_pathogen_p)>=jsd_th & diff_ctrl_test>tol_in)
  ids_matter    = unique(ids_matter_df %>% dplyr::pull(param_set_id))
  length(ids_matter)
  df_comparisons = df_comparisons %>% dplyr::filter(param_set_id %in% ids_matter)
}

# ============= MACSPEC-SPECIFIC ANALYSIS =============
if(perform_macspec_analysis) {
  df_results_pick_p = df_results %>%
    dplyr::filter(injury_type=='pathogenic' &
                    get(paste0('d_macspec_vs_', group2, '_pathogen_e'))>0.3)
  df_results_pick_s = df_results %>%
    dplyr::filter(injury_type=='sterile' &
                    get(paste0('d_macspec_vs_', group2, '_sterile_e'))>0.3)

  print(paste("Mean", group1, "pathogenic:",
              mean(df_results_pick_p[[paste0('mean_', group1, '_pathogen_e')]])))
  print(paste("Mean", group2, "pathogenic:",
              mean(df_results_pick_p[[paste0('mean_', group2, '_pathogen_e')]])))
  print(paste("Mean", group1, "sterile:",
              mean(df_results_pick_s[[paste0('mean_', group1, '_sterile_e')]])))
  print(paste("Mean", group2, "sterile:",
              mean(df_results_pick_s[[paste0('mean_', group2, '_sterile_e')]])))
}

# ============= CALCULATE DIFFERENCES =============
# Construct column names dynamically
group1_pathogen_col = paste0('mean_', group1, '_pathogen_e')
group1_sterile_col  = paste0('mean_', group1, '_sterile_e')
group2_pathogen_col = paste0('mean_', group2, '_pathogen_e')
group2_sterile_col  = paste0('mean_', group2, '_sterile_e')

df_comparisons = df_comparisons %>%
  dplyr::mutate(diff_group1_minus_group2 = ifelse(injury_type=='pathogenic',
                                                   .data[[group1_pathogen_col]] - .data[[group2_pathogen_col]],
                                                   .data[[group1_sterile_col]] - .data[[group2_sterile_col]]))

# ============= RANDOMNESS CHECK (if requested) =============
if(perform_randomness_check) {
  df_comparisons = df_comparisons %>%
    dplyr::mutate(diff_notrnd_minus_rnd = ifelse(injury_type=='pathogenic',
                                                  mean_tregs_on_pathogen_e - mean_tregs_rnd_pathogen_e,
                                                  mean_tregs_on_sterile_e - mean_tregs_rnd_sterile_e))
}

# ============= MERGE WITH COHEN'S D =============
# Construct Cohen's d column names
cohens_d_pathogen_col = paste0('d_', group1, '_vs_', group2, '_pathogen_e')
cohens_d_sterile_col  = paste0('d_', group1, '_vs_', group2, '_sterile_e')

# Handle alternative naming conventions for Cohen's d columns
# For macspec comparisons, the column names might be reversed
if(!cohens_d_pathogen_col %in% colnames(df_results)) {
  cohens_d_pathogen_col_alt = paste0('d_', group2, '_vs_', group1, '_pathogen_e')
  cohens_d_sterile_col_alt  = paste0('d_', group2, '_vs_', group1, '_sterile_e')
  if(cohens_d_pathogen_col_alt %in% colnames(df_results)) {
    cohens_d_pathogen_col = cohens_d_pathogen_col_alt
    cohens_d_sterile_col  = cohens_d_sterile_col_alt
  }
}

select_cols = c('param_set_id', 'diff_group1_minus_group2')
merge_cols  = c('param_set_id', cohens_d_pathogen_col, cohens_d_sterile_col)

if(perform_randomness_check) {
  select_cols = c(select_cols, 'diff_notrnd_minus_rnd')
  merge_cols  = c(merge_cols, 'd_tregs_rnd_vs_nrnd_pathogen_e', 'd_tregs_rnd_vs_nrnd_sterile_e')
}

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, all_of(select_cols[-1]))
df_comparisons = merge(df_comparisons,
                       distinct(df_results[merge_cols]),
                       by='param_set_id')

# ============= AUTO-CALCULATE TOLERANCE IF NEEDED =============
if(is.na(tol_in)) {
  df_comparisons_below = df_comparisons %>%
    dplyr::filter((injury_type=='pathogenic' & abs(.data[[cohens_d_pathogen_col]])<jsd_th) |
                  (injury_type=='sterile' & abs(.data[[cohens_d_sterile_col]])<jsd_th))
  tol_in = max(max(df_comparisons_below$diff_group1_minus_group2),
               abs(min(df_comparisons_below$diff_group1_minus_group2)))
  print(paste("Auto-calculated tol_in:", tol_in))
}

# ============= RANDOMNESS SENSE CHECK =============
if(perform_randomness_check) {
  df_comparisons_keep  = df_comparisons
  df_comparisons_sense = df_comparisons %>%
    dplyr::mutate(makes_sense = ifelse(
      sign(diff_group1_minus_group2) != sign(diff_notrnd_minus_rnd) &
        abs(diff_notrnd_minus_rnd) > tol_in &
        (abs(d_tregs_rnd_vs_nrnd_pathogen_e) >= jsd_th | abs(d_tregs_rnd_vs_nrnd_sterile_e) >= jsd_th),
      0, 1))
  df_comparisons_sense = df_comparisons_sense %>% dplyr::filter(makes_sense==0)
  print(paste("Number of cases that don't make sense:", dim(df_comparisons_sense)[1]))
  df_comparisons = df_comparisons_keep
}

# ============= PREPARE PLOT DATA =============
df_comparisons_plot = df_comparisons

# --- PATHOGENIC ---
df_comparisons_plot_pathogenic = df_comparisons_plot %>%
  dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>%
  dplyr::mutate(tregs_better = ifelse(diff_group1_minus_group2 > tol_in, 1,
                                      ifelse(diff_group1_minus_group2 < -1*tol_in, -1, 0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>%
  dplyr::mutate(tregs_better_cohens = ifelse(abs(.data[[cohens_d_pathogen_col]]) > jsd_th,
                                             tregs_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>%
  dplyr::select(-all_of(cohens_d_sterile_col))
colnames(df_comparisons_plot_pathogenic)[which(
  colnames(df_comparisons_plot_pathogenic)==cohens_d_pathogen_col)] = 'cohens_d'

# --- STERILE ---
df_comparisons_plot_sterile = df_comparisons_plot %>%
  dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>%
  dplyr::mutate(tregs_better = ifelse(diff_group1_minus_group2 > tol_in, 1,
                                      ifelse(diff_group1_minus_group2 < -1*tol_in, -1, 0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>%
  dplyr::mutate(tregs_better_cohens = ifelse(abs(.data[[cohens_d_sterile_col]]) > jsd_th,
                                             tregs_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>%
  dplyr::select(-all_of(cohens_d_pathogen_col))
colnames(df_comparisons_plot_sterile)[which(
  colnames(df_comparisons_plot_sterile)==cohens_d_sterile_col)] = 'cohens_d'

# ============= CONFLICT ANALYSIS =============
tregs_better_sterile    = df_comparisons_plot_sterile %>%
  dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_sterile     = df_comparisons_plot_sterile %>%
  dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)
tregs_better_pathogenic = df_comparisons_plot_pathogenic %>%
  dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_pathogenic  = df_comparisons_plot_pathogenic %>%
  dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)

s_better_p_better = intersect(tregs_better_sterile, tregs_better_pathogenic)
s_better_p_worse  = intersect(tregs_better_sterile, tregs_worse_pathogenic)
s_worse_p_better  = intersect(tregs_worse_sterile, tregs_better_pathogenic)
s_worse_p_worse   = intersect(tregs_worse_sterile, tregs_worse_pathogenic)

# ============= COMBINE AND SAVE =============
df_comparisons_plot = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)

# Rename the difference column for plotting
colnames(df_comparisons_plot)[which(
  colnames(df_comparisons_plot)=='diff_group1_minus_group2')] = 'diff_tregs_on_minus_off'

# Save results
output_rds = paste0('./df_comparisons_plot_', comparison_name, '.rds')
saveRDS(df_comparisons_plot, output_rds)
print(paste("Saved results to:", output_rds))

# ============= GENERATE PLOT =============
cohens_th = jsd_th
e_th      = tol_in

source('./MISC/REGIONS.R')

output_plot = paste0("./ABM_JSD_", toupper(comparison_name), ".png")
ggsave(
  filename = output_plot,
  plot = p,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)
print(paste("Saved plot to:", output_plot))

# ============= SUMMARY STATISTICS =============
print(paste("Total unique parameter sets:", length(unique(dfp$param_set_id))))

dfp_p = dfp %>% dplyr::filter(injury_type=='pathogenic') %>%
  dplyr::select(param_set_id, tregs_better_cohens)
dfp_s = dfp %>% dplyr::filter(injury_type=='sterile') %>%
  dplyr::select(param_set_id, tregs_better_cohens)
colnames(dfp_p)[2] = 'tregs_better_cohens_p'
colnames(dfp_s)[2] = 'tregs_better_cohens_s'
dfp_sp = merge(dfp_p, dfp_s, by='param_set_id')
dfp_sp = dfp_sp %>%
  dplyr::mutate(treg_better_cohens_both = ifelse(tregs_better_cohens_p + tregs_better_cohens_s == 2, 1, 0))

print(paste("% better for pathogenic:",
            round(100*length(which(dfp_sp$tregs_better_cohens_p==1))/dim(dfp_sp)[1], 2)))
print(paste("% worse for pathogenic:",
            round(100*length(which(dfp_sp$tregs_better_cohens_p==-1))/dim(dfp_sp)[1], 2)))
print(paste("% better for sterile:",
            round(100*length(which(dfp_sp$tregs_better_cohens_s==1))/dim(dfp_sp)[1], 2)))
print(paste("% worse for sterile:",
            round(100*length(which(dfp_sp$tregs_better_cohens_s==-1))/dim(dfp_sp)[1], 2)))
print(paste("% better for both:",
            round(100*sum(dfp_sp$treg_better_cohens_both)/dim(dfp_sp)[1], 2)))

df_comparisons_plot_nz = df_comparisons_plot %>% dplyr::filter(tregs_better_cohens!=0)

# Check conflicts
print("Conflicting cases (better in one, worse in other):")
print(c(s_better_p_worse, s_worse_p_better))
