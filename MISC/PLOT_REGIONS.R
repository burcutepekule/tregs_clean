if(score_type=='epithelial'){
  e_th = tol_in_e
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
}else if(score_type=='pathogen'){
  e_th = tol_in_p
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
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > e_th, -1,
                                                                                                         ifelse(diff_compare < -1*e_th,+1,0)))
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                                diff_better, 0))
  
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
  colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'
  
  # ====================== STERILE ==================================================================
  df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
  df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > e_th, -1,
                                                                                                   ifelse(diff_compare < -1*e_th,+1,0)))
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
  df_comparisons_plot$diff_compare = -1*df_comparisons_plot$diff_compare  # for plotting purposes
}else{
  df_comparisons_keep = df_comparisons
  # ==========================================================================================
  score_type_use = 'epithelial'
  mean_subt_from_full_p = paste0('mean_',condition_subt_from,'_pathogen_',substr(score_type_use,1,1))
  mean_subt_from_full_s = paste0('mean_',condition_subt_from,'_sterile_',substr(score_type_use,1,1))
  
  mean_subt_full_p = paste0('mean_',condition_subt,'_pathogen_',substr(score_type_use,1,1))
  mean_subt_full_s = paste0('mean_',condition_subt,'_sterile_',substr(score_type_use,1,1))
  
  d_full_p = paste0('d_',jensen_distance,'_pathogen_',substr(score_type_use,1,1))
  d_full_s = paste0('d_',jensen_distance,'_sterile_',substr(score_type_use,1,1))
  
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
  
  df_comparisons_plot_pathogenic$score_type_use = score_type_use
  df_comparisons_plot_sterile$score_type_use    = score_type_use
  df_comparisons_plot_1 = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)
  
  # ==========================================================================================
  df_comparisons = df_comparisons_keep
  score_type_use = 'pathogen'
  # ==========================================================================================
  mean_subt_from_full_p = paste0('mean_',condition_subt_from,'_pathogen_',substr(score_type_use,1,1))
  mean_subt_from_full_s = paste0('mean_',condition_subt_from,'_sterile_',substr(score_type_use,1,1))
  
  mean_subt_full_p = paste0('mean_',condition_subt,'_pathogen_',substr(score_type_use,1,1))
  mean_subt_full_s = paste0('mean_',condition_subt,'_sterile_',substr(score_type_use,1,1))
  
  d_full_p = paste0('d_',jensen_distance,'_pathogen_',substr(score_type_use,1,1))
  d_full_s = paste0('d_',jensen_distance,'_sterile_',substr(score_type_use,1,1))
  
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
  # ====================== PATHOGENIC ==================================================================
  df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_p, -1,
                                                                                                         ifelse(diff_compare < -1*tol_in_p,+1,0)))
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_p]])>jsd_th,
                                                                                                                diff_better, 0))
  
  df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-all_of(d_full_s))
  colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)==d_full_p)]='cohens_d'
  
  # ====================== STERILE ==================================================================
  df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
  df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better = ifelse(diff_compare > tol_in_p, -1,
                                                                                                   ifelse(diff_compare < -1*tol_in_p,+1,0)))
  df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(diff_better_cohens = ifelse(abs(.data[[d_full_s]])>jsd_th,
                                                                                                          diff_better, 0))
  
  df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-all_of(d_full_p))
  colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)==d_full_s)]='cohens_d'
  
  diff_better_sterile_p_score    = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
  diff_worse_sterile_p_score     = df_comparisons_plot_sterile %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
  diff_better_pathogenic_p_score = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
  diff_worse_pathogenic_p_score  = df_comparisons_plot_pathogenic %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
  
  # ============CONFLICT?=========================================================
  df_comparisons_plot_pathogenic$score_type_use = score_type_use
  df_comparisons_plot_sterile$score_type_use    = score_type_use
  df_comparisons_plot_2 = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)
  
  df_comparisons_plot_12 = rbind(df_comparisons_plot_1, df_comparisons_plot_2)
  
  df_combined = df_comparisons_plot_12 %>%
    # Pivot wider to get epithelial and pathogenic side by side
    pivot_wider(
      id_cols = c(param_set_id, injury_type),
      names_from = score_type_use,
      values_from = c(diff_compare, cohens_d, diff_better, diff_better_cohens)
    ) 
  
  df_comparisons_plot = df_combined %>% dplyr::rowwise() %>%
    dplyr::mutate(
      # New diff_better_cohens -> categorize for higher e score + lower p score (better in both cases)
      # or the opposite (worse for both cases)
      diff_better_cohens = case_when(
        diff_better_cohens_epithelial == 1 & diff_better_cohens_pathogen == 1 ~ 1,
        diff_better_cohens_epithelial == -1 & diff_better_cohens_pathogen == -1 ~ -1,
        TRUE ~ 0
      )
    )
  
  df_comparisons_plot = df_comparisons_plot %>% dplyr::rowwise() %>% dplyr::mutate(
    # Sum of Cohen's d values
    # cohens_d = 0.5*(cohens_d_epithelial+cohens_d_pathogen),
    cohens_d = min(cohens_d_epithelial, cohens_d_pathogen)
  )

  df_comparisons_plot = df_comparisons_plot %>% dplyr::rowwise() %>% dplyr::mutate(

      diff_compare_e = case_when(
        diff_compare_epithelial > tol_in_e ~ diff_compare_epithelial-tol_in_e,
        diff_compare_epithelial < -1*tol_in_e ~ diff_compare_epithelial+tol_in_e,
        TRUE ~ 0
      ),
      
      diff_compare_p = case_when(
        diff_compare_pathogen > tol_in_p ~ (-1*diff_compare_pathogen+tol_in_p)*0.01,
        diff_compare_pathogen < -1*tol_in_p ~ (-1*diff_compare_pathogen-tol_in_p)*0.01,
        TRUE ~ 0
      ),
      
      diff_compare = if_else(
        abs(diff_compare_e) < abs(diff_compare_p),
        diff_compare_e,
        diff_compare_p
      )
    ) 
  
  df_comparisons_plot = df_comparisons_plot[c("param_set_id","injury_type","diff_compare","cohens_d","diff_better_cohens")]
  
  ### FILTER FOR PATHOGENIC BECAUSE "BOTH" MEANS THERE NEEDS TO BE PATHOGEN ABUNDANCE AS WELL!
  df_comparisons_plot = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
  # e_th                = tol_in_p # example threshold for y
  e_th                = 0 # example threshold for y
  
}

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
  }else if(num_cols>12){
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
  filename = paste0("./ABM_JSD_",jensen_distance,"_",score_type,".png"),
  plot = p_use,
  width = width_adjust,
  height = 6,
  dpi = 300,
  bg='white'
)
