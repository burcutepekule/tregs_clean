# is_under_control = function(epithelial_health, pathogen_load) {
#   return(epithelial_health > 150*0.75 & pathogen_load < 10)
# }

is_under_control = function(epithelial_health, pathogen_load, epithelial_limit, pathogen_limit) {
  return(epithelial_health < epithelial_limit & pathogen_load < pathogen_limit)
}

# Clean up dynamically created variables from previous iterations
rm(list = ls(pattern = "^(scores_|time_ss_|full_data_comparison_scores_|osc_)"))

# message(paste0("Re-processing param_set_", param_id,", optimization index ",i_opt,", for overwrite_in=", overwrite_in))

results_list = lapply(ros_vals, function(ros_in) {
  lapply(pat_vals, function(pat_in) {
    file_path_data = paste0(path_data, 'longitudinal_df_param_set_id_', param_id, 
                            '_sterile_', sterile_in,
                            '_macspec_', macspec_in,
                            '_tregs_', tregs_on_in,
                            '_ros_level_', ros_in,
                            '_pat_level_', pat_in,
                            '_trnd_', trnd_in,
                            '_overwrite_', overwrite_in,
                            '_optidx_', i_opt,
                            '_effidx_', treg_eff_in,
                            '_m2on_', m2on_in,
                            '.rds')
    
    if(file.exists(file_path_data)) {
      return(readRDS(file_path_data))
    }
    return(NULL)
  })
})

if (is.null(results_list[[1]][[1]])) {
  message(paste0("Skipping ", param_id))
  skip_this = TRUE
}else {
  skip_this = FALSE
  # Flatten and combine in one step
  results = rbindlist(unlist(results_list, recursive = FALSE), use.names = TRUE, fill = TRUE)
  results = as.data.frame(results)
  
  results_all = results %>%
    dplyr::group_by(param_set_id, sterile, macspec_on, tregs_on,
                    randomize_tregs, ros_level, pat_level, rep_id) %>%
    dplyr::summarise(mean_score_e_all = mean(epithelial_score, na.rm = TRUE), 
                     mean_score_p_all = mean(pathogen, na.rm = TRUE), .groups = "drop")
  
  
  results_150 = results %>% dplyr::filter(t>=(max(results$t)-150))
  
  results_150 = results_150 %>%
    dplyr::group_by(param_set_id, sterile, macspec_on, tregs_on,
                    randomize_tregs, ros_level, pat_level, rep_id) %>%
    dplyr::summarise(mean_score_e = mean(epithelial_score, na.rm = TRUE), 
                     mean_score_p = mean(pathogen, na.rm = TRUE), .groups = "drop")
  
  results_150 = merge(results_150, results_all)
  
  results_150 = results_150 %>% dplyr::filter(!(mean_score_p_all==0 & mean_score_e_all==0))
  
  
  results_150 = results_150 %>%
    pivot_longer(
      cols = c(mean_score_e, mean_score_p),
      names_to = "score_type",
      values_to = "mean_score",
      names_prefix = "mean_score_"
    ) %>%
    mutate(score_type = recode(score_type, "e" = "epithelium", "p" = "pathogen")) %>% rename(replicate_id = rep_id)
  
  print('done filtering 150')
  
  
  comparison_results = results_150
  comparison_results$injury_type = 'pathogenic'
  comparison_results$macspec_on  = 0
  comparison_results$tregs_on    = tregs_on_in
  comparison_results$tregs_rnd   = 0
  
  all_comparison_results = comparison_results
  
  df_steps = all_comparison_results %>% dplyr::filter(injury_type == 'pathogenic')
  
  # Create separate dataframes for epithelium and pathogen scores
  df_epithelium = df_steps %>%
    dplyr::filter(score_type == "epithelium") %>%
    dplyr::select(param_set_id, replicate_id, ros_level, pat_level, mean_score) %>%
    dplyr::rename(mean_epithelium = mean_score)
  
  df_pathogen = df_steps %>%
    dplyr::filter(score_type == "pathogen") %>%
    dplyr::select(param_set_id, replicate_id, ros_level, pat_level, mean_score) %>%
    dplyr::rename(mean_pathogen = mean_score)
  
  # Join epithelium and pathogen data for each replicate
  df_combined = df_epithelium %>%
    dplyr::inner_join(df_pathogen,
                      by = c("param_set_id", "replicate_id", "ros_level", "pat_level"))
  
  # Calculate per-replicate control status using the per-replicate mean scores
  df_combined = df_combined %>%
    dplyr::mutate(is_controlled = is_under_control(mean_epithelium, mean_pathogen,
                                                   epithelial_limit, pathogen_limit))
  
  print('computing control matrix')
  
  # Build the control matrix (now as percentage of controlled replicates)
  control_matrix = matrix(NA, nrow = length(pat_vals), ncol = length(ros_vals))
  rownames(control_matrix) = paste0("pat", pat_vals)
  colnames(control_matrix) = paste0("ros", ros_vals)
  
  control_matrix_long = c()
  for (i in seq_along(pat_vals)) {
    for (j in seq_along(ros_vals)) {
      pat = pat_vals[i]
      ros = ros_vals[j]
      
      df_combined %>%
        dplyr::filter(param_set_id == param_id,
                      ros_level == ros,
                      pat_level == pat)
      
      # Get all control status values for this ros/pat combination across replicates
      controlled_vals = df_combined %>%
        dplyr::filter(param_set_id == param_id,
                      ros_level == ros,
                      pat_level == pat) %>%
        pull(is_controlled)
      
      # Calculate percentage of replicates that are controlled
      if(length(controlled_vals) > 0) {
        control_matrix[i, j] = sum(controlled_vals) / length(controlled_vals)
        control_matrix_long  = rbind(control_matrix_long, c(param_id, pat, ros, sum(controlled_vals) / length(controlled_vals)))
      } else {
        control_matrix[i, j] = 0
      }
    }
  }
}