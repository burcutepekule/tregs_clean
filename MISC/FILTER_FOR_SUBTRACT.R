e_th = tol_in_e
# ==========================================================================================
mean_subt_from_full_p = paste0('mean_',condition_subt_from,'_pathogen_',substr(score_type,1,1))
mean_subt_from_full_s = paste0('mean_',condition_subt_from,'_sterile_',substr(score_type,1,1))

mean_subt_full_p = paste0('mean_',condition_subt,'_pathogen_',substr(score_type,1,1))
mean_subt_full_s = paste0('mean_',condition_subt,'_sterile_',substr(score_type,1,1))

d_full_p = paste0('d_',jensen_distance,'_pathogen_',substr(score_type,1,1))
d_full_s = paste0('d_',jensen_distance,'_sterile_',substr(score_type,1,1))

df_comparisons = df_comparisons_in %>% 
  dplyr::mutate(diff_compare = ifelse(
    injury_type == 'pathogenic', 
    .data[[mean_subt_from_full_p]] - .data[[mean_subt_full_p]],
    .data[[mean_subt_from_full_s]] - .data[[mean_subt_full_s]]
  ))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_compare)
df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id',d_full_p, d_full_s)]), by='param_set_id')

# ====================== Find significant cases ======================================================
df_comparisons_plot            = df_comparisons
# ====================== PATHOGENIC ==================================================================
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > e_th, 1,
                                                                                                       ifelse(diff_compare < -1*e_th,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                              diff_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'

# ====================== STERILE ==================================================================
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > e_th, 1,
                                                                                                 ifelse(diff_compare < -1*e_th,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_s]])>jsd_th,
                                                                                                        diff_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-all_of(d_full_p))
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)==d_full_s)]='cohens_d'

df_comparisons_plot = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)