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

# ============CONFLICT?=========================================================
diff_better_sterile    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_sterile     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
diff_better_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_pathogenic  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)

s_better_p_better = intersect(diff_better_sterile, diff_better_pathogenic)
s_better_p_worse  = intersect(diff_better_sterile, tregs_worse_pathogenic)
s_worse_p_better  = intersect(tregs_worse_sterile, diff_better_pathogenic)
s_worse_p_worse   = intersect(tregs_worse_sterile, tregs_worse_pathogenic)
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

cohens_th           = jsd_th # example threshold for x
e_th                = tol_in # example threshold for y

source('./MISC/REGIONS.R')

if(exists("num_cols")){
  if(num_cols<5){
    width_adjust = 9
    if(labels_on==1){
      p_use = p_label_on
    }else{
      p_use = p_label_off
    }
  }else if(num_cols>12){
    width_adjust = 21
    labels_on    = 0 # flip because too many
    p_use = p_label_off
  }else{
    width_adjust = round(num_cols*1.75)
    if(labels_on==1){
      p_use = p_label_on
    }else{
      p_use = p_label_off
    }
  }
}else{
  width_adjust = 9
  p_use = p_label_off
}


ggsave(
  filename = paste0("./ABM_JSD_",jensen_distance,"_",score_type,".png"),
  plot = p_use,
  width = width_adjust,
  height = 6,
  dpi = 300,
  bg='white'
)
