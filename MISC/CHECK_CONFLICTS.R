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
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_e, 1,
                                                                                                       ifelse(diff_compare < -1*tol_in_e,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                              diff_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'

# ====================== STERILE ==================================================================
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_e, 1,
                                                                                                 ifelse(diff_compare < -1*tol_in_e,-1,0)))
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
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_e, 1,
                                                                                                       ifelse(diff_compare < -1*tol_in_e,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                              diff_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'

# ====================== STERILE ==================================================================
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_e, 1,
                                                                                                 ifelse(diff_compare < -1*tol_in_e,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_s]])>jsd_th,
                                                                                                        diff_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-all_of(d_full_p))
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)==d_full_s)]='cohens_d'

# HERE -1 +1 SHOULD SWITCH BECAUSE HIGHER PATHOGEN MEANS WORSE FOR PATHOGEN CONTOL! IT'S NOT LIKE EPITHELIAL SCORE! 
diff_better_sterile_p_score    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_worse_sterile_p_score     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==+1) %>% dplyr::pull(param_set_id)
diff_better_pathogenic_p_score = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_worse_pathogenic_p_score  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==+1) %>% dplyr::pull(param_set_id)

