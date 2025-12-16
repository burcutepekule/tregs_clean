if(score_type=='epithelial'){
  e_th = tol_in_e
}else if(score_type=='pathogen'){
  e_th = tol_in_p
}

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

# ============CONFLICT?=========================================================
diff_better_sterile    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
diff_worse_sterile     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_better_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
diff_worse_pathogenic  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)

s_better_p_better = intersect(diff_better_sterile, diff_better_pathogenic)
s_better_p_worse  = intersect(diff_better_sterile, diff_worse_pathogenic)
s_worse_p_better  = intersect(diff_worse_sterile, diff_better_pathogenic)
s_worse_p_worse   = intersect(diff_worse_sterile, diff_worse_pathogenic)
conflicting       = c(s_better_p_worse, s_worse_p_better)
# ==============================================================================

df_comparisons_plot = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)

if(length(conflicting)>0){
  df_comparisons_plot_conflicting     = df_comparisons_plot %>% dplyr::filter(param_set_id %in% conflicting)
  df_comparisons_plot_not_conflicting = df_comparisons_plot %>% dplyr::filter(!(param_set_id %in% conflicting))
  df_comparisons_plot_conflicting     = df_comparisons_plot_conflicting %>% 
    dplyr::mutate(param_set_id = paste0(param_set_id,'_',substr(injury_type, 1, 1)))
  df_comparisons_plot = rbind(df_comparisons_plot_conflicting, df_comparisons_plot_not_conflicting)
}
df_comparisons_plot = df_comparisons_plot[c("param_set_id","injury_type","diff_compare","cohens_d","diff_better_cohens")]

cohens_th = jsd_th # example threshold for x

source('./MISC/REGIONS.R')

if(exists("num_cols")){
  if(num_cols<5){
    width_adjust = 9
    if(labels_on==1){
      p_use = p_label_on
    }else{
      p_use = p_label_off
    }
  }else if(num_cols>18){
    width_adjust = 9
    labels_on    = 0 # flip because too many
    p_use = p_label_off
  }else{
    if(labels_on==1){
      width_adjust = round(num_cols*1.75)
      p_use = p_label_on
    }else{
      width_adjust = 9
      p_use = p_label_off
    }
  }
}else{
  width_adjust = 9
  p_use = p_label_off
}

ggsave(
  filename = paste0("./ABM_JSD_",jensen_distance,"_",score_type,"_",filter_control,".png"),
  plot = p_use,
  width = width_adjust,
  height = 6,
  dpi = 300,
  bg='white'
)
