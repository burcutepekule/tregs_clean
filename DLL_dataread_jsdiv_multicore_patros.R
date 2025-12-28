rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(philentropy)

# Get the directory of the current script
get_script_dir = function() {
  # Try rstudioapi first (works in RStudio)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  # Try commandArgs (works with Rscript)
  args = commandArgs(trailingOnly = FALSE)
  file_arg = grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # Try sys.frame (works with source())
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(sys.frame(1)$ofile))
  }
  # Fallback to current working directory
  return(getwd())
}

setwd(get_script_dir())
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Load files
results_merged             = c()
sterile_comparison_keep    = c()
pathogenic_comparison_keep = c()

path = "/Users/burcutepekule/Desktop/sim_abm/"
# File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
files_0_1 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds$", full.names = TRUE)
files_0_2 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_2_trnd_0.rds$", full.names = TRUE)
files_0_3 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_3_trnd_0.rds$", full.names = TRUE)
files_1_1 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_1_trnd_0.rds$", full.names = TRUE)
files_1_2 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_2_trnd_0.rds$", full.names = TRUE)
files_1_3 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_3_trnd_0.rds$", full.names = TRUE)
files_2_1 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_1_trnd_0.rds$", full.names = TRUE)
files_2_2 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_2_trnd_0.rds$", full.names = TRUE)
files_2_3 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_3_trnd_0.rds$", full.names = TRUE)

indices_0_1 = str_extract(basename(files_0_1), "\\d+") |> as.numeric()
indices_0_2 = str_extract(basename(files_0_2), "\\d+") |> as.numeric()
indices_0_3 = str_extract(basename(files_0_3), "\\d+") |> as.numeric()
indices_1_1 = str_extract(basename(files_1_1), "\\d+") |> as.numeric()
indices_1_2 = str_extract(basename(files_1_2), "\\d+") |> as.numeric()
indices_1_3 = str_extract(basename(files_1_3), "\\d+") |> as.numeric()
indices_2_1 = str_extract(basename(files_2_1), "\\d+") |> as.numeric()
indices_2_2 = str_extract(basename(files_2_2), "\\d+") |> as.numeric()
indices_2_3 = str_extract(basename(files_2_3), "\\d+") |> as.numeric()

indices = Reduce(intersect, list(
  indices_0_1, indices_0_2, indices_0_3,
  indices_1_1, indices_1_2, indices_1_3,
  indices_2_1, indices_2_2, indices_2_3
))

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('./ids_read_abm_patros.rds')){
  inds_read = readRDS('./ids_read_abm_patros.rds')
}else{
  inds_read = c()
}

inds2read = sort(setdiff(indices,inds_read))
length(inds2read)

split_equal = function(x, n_chunks) {
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

# n1 = 1
# n2 = 1

if(n1==1){
  loop_over = inds2read
}else{
  chunks    = split_equal(inds2read, n1)
  loop_over = chunks[[n2]]
}

if(length(loop_over)>0){
  # Track indices that have been successfully processed
  processed_indices = c()
  
  for (i_idx in seq_along(loop_over)){
    
    i = loop_over[i_idx]
    message("Processing param_set_", i)
    
    files_to_check = c(
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_2_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_3_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_2_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_3_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_2_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_3_trnd_0.rds')
    )
    if(any(file.info(files_to_check)$size<100000)){
      processed_indices      = c(processed_indices, i) #add and skip
      message("Skipped one")
    }else{
      
      results_0_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds'))
      results_0_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_2_trnd_0.rds'))
      results_0_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_3_trnd_0.rds'))
      results_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_1_trnd_0.rds'))
      results_1_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_2_trnd_0.rds'))
      results_1_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_1_pat_level_3_trnd_0.rds'))
      results_2_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_1_trnd_0.rds'))
      results_2_2 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_2_trnd_0.rds'))
      results_2_3 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_2_pat_level_3_trnd_0.rds'))
      
      results = rbind(
        results_0_1, results_0_2, results_0_3,
        results_1_1, results_1_2, results_1_3,
        results_2_1, results_2_2, results_2_3
      )
      
      full_data_comparison = results %>% dplyr::select(param_set_id, sterile, macspec_on, tregs_on, 
                                                       randomize_tregs, ros_level, pat_level, rep_id, 
                                                       t, time_ss_e, time_ss_p, epithelial_score, pathogen)
      min_reps  = min(full_data_comparison$rep_id)
      max_reps  = max(full_data_comparison$rep_id)
      t_max_ind = max(full_data_comparison$t)
      
      # ====== PATHOGEN
      # Test scenarios: sterile_0_macspec_0_tregs_0_ros_level_{0,1,2}_pat_level_{1,2,3}_trnd_0
      scores_0_1_p_keep = c()
      scores_0_2_p_keep = c()
      scores_0_3_p_keep = c()
      scores_1_1_p_keep = c()
      scores_1_2_p_keep = c()
      scores_1_3_p_keep = c()
      scores_2_1_p_keep = c()
      scores_2_2_p_keep = c()
      scores_2_3_p_keep = c()
      
      # === EPITHELIUM
      # Test scenarios: sterile_0_macspec_0_tregs_0_ros_level_{0,1,2}_pat_level_{1,2,3}_trnd_0
      scores_0_1_e_keep = c()
      scores_0_2_e_keep = c()
      scores_0_3_e_keep = c()
      scores_1_1_e_keep = c()
      scores_1_2_e_keep = c()
      scores_1_3_e_keep = c()
      scores_2_1_e_keep = c()
      scores_2_2_e_keep = c()
      scores_2_3_e_keep = c()
      
      all_comparison_results_reps = data.frame()
      
      for (rep in min_reps:max_reps) {
        
        full_data_comparison_scores_0_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==0 & pat_level==1 & randomize_tregs==0)
        full_data_comparison_scores_0_2 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==0 & pat_level==2 & randomize_tregs==0)
        full_data_comparison_scores_0_3 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==0 & pat_level==3 & randomize_tregs==0)
        full_data_comparison_scores_1_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==1 & pat_level==1 & randomize_tregs==0)
        full_data_comparison_scores_1_2 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==1 & pat_level==2 & randomize_tregs==0)
        full_data_comparison_scores_1_3 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==1 & pat_level==3 & randomize_tregs==0)
        full_data_comparison_scores_2_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==2 & pat_level==1 & randomize_tregs==0)
        full_data_comparison_scores_2_2 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==2 & pat_level==2 & randomize_tregs==0)
        full_data_comparison_scores_2_3 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 & ros_level==2 & pat_level==3 & randomize_tregs==0)
        
        time_ss_0_1_e = as.numeric(steady_state_idx(full_data_comparison_scores_0_1$epithelial_score))
        time_ss_0_2_e = as.numeric(steady_state_idx(full_data_comparison_scores_0_2$epithelial_score))
        time_ss_0_3_e = as.numeric(steady_state_idx(full_data_comparison_scores_0_3$epithelial_score))
        time_ss_1_1_e = as.numeric(steady_state_idx(full_data_comparison_scores_1_1$epithelial_score))
        time_ss_1_2_e = as.numeric(steady_state_idx(full_data_comparison_scores_1_2$epithelial_score))
        time_ss_1_3_e = as.numeric(steady_state_idx(full_data_comparison_scores_1_3$epithelial_score))
        time_ss_2_1_e = as.numeric(steady_state_idx(full_data_comparison_scores_2_1$epithelial_score))
        time_ss_2_2_e = as.numeric(steady_state_idx(full_data_comparison_scores_2_2$epithelial_score))
        time_ss_2_3_e = as.numeric(steady_state_idx(full_data_comparison_scores_2_3$epithelial_score))
        
        time_ss_0_1_p = as.numeric(steady_state_idx(full_data_comparison_scores_0_1$pathogen))
        time_ss_0_2_p = as.numeric(steady_state_idx(full_data_comparison_scores_0_2$pathogen))
        time_ss_0_3_p = as.numeric(steady_state_idx(full_data_comparison_scores_0_3$pathogen))
        time_ss_1_1_p = as.numeric(steady_state_idx(full_data_comparison_scores_1_1$pathogen))
        time_ss_1_2_p = as.numeric(steady_state_idx(full_data_comparison_scores_1_2$pathogen))
        time_ss_1_3_p = as.numeric(steady_state_idx(full_data_comparison_scores_1_3$pathogen))
        time_ss_2_1_p = as.numeric(steady_state_idx(full_data_comparison_scores_2_1$pathogen))
        time_ss_2_2_p = as.numeric(steady_state_idx(full_data_comparison_scores_2_2$pathogen))
        time_ss_2_3_p = as.numeric(steady_state_idx(full_data_comparison_scores_2_3$pathogen))
        
        time_ss_vec = c(
          time_ss_0_1_e, time_ss_0_1_p,
          time_ss_0_2_e, time_ss_0_2_p,
          time_ss_0_3_e, time_ss_0_3_p,
          time_ss_1_1_e, time_ss_1_1_p,
          time_ss_1_2_e, time_ss_1_2_p,
          time_ss_1_3_e, time_ss_1_3_p,
          time_ss_2_1_e, time_ss_2_1_p,
          time_ss_2_2_e, time_ss_2_2_p,
          time_ss_2_3_e, time_ss_2_3_p
        )
        
        if(!any(is.na(time_ss_vec))){
          
          # ==== PATHOGEN ABUNDANCE
          scores_0_1_p = full_data_comparison_scores_0_1$pathogen[time_ss_0_1_p:t_max_ind]
          scores_0_2_p = full_data_comparison_scores_0_2$pathogen[time_ss_0_2_p:t_max_ind]
          scores_0_3_p = full_data_comparison_scores_0_3$pathogen[time_ss_0_3_p:t_max_ind]
          scores_1_1_p = full_data_comparison_scores_1_1$pathogen[time_ss_1_1_p:t_max_ind]
          scores_1_2_p = full_data_comparison_scores_1_2$pathogen[time_ss_1_2_p:t_max_ind]
          scores_1_3_p = full_data_comparison_scores_1_3$pathogen[time_ss_1_3_p:t_max_ind]
          scores_2_1_p = full_data_comparison_scores_2_1$pathogen[time_ss_2_1_p:t_max_ind]
          scores_2_2_p = full_data_comparison_scores_2_2$pathogen[time_ss_2_2_p:t_max_ind]
          scores_2_3_p = full_data_comparison_scores_2_3$pathogen[time_ss_2_3_p:t_max_ind]
          
          # Accumulate scores
          scores_0_1_p_keep = c(scores_0_1_p_keep, scores_0_1_p)
          scores_0_2_p_keep = c(scores_0_2_p_keep, scores_0_2_p)
          scores_0_3_p_keep = c(scores_0_3_p_keep, scores_0_3_p)
          scores_1_1_p_keep = c(scores_1_1_p_keep, scores_1_1_p)
          scores_1_2_p_keep = c(scores_1_2_p_keep, scores_1_2_p)
          scores_1_3_p_keep = c(scores_1_3_p_keep, scores_1_3_p)
          scores_2_1_p_keep = c(scores_2_1_p_keep, scores_2_1_p)
          scores_2_2_p_keep = c(scores_2_2_p_keep, scores_2_2_p)
          scores_2_3_p_keep = c(scores_2_3_p_keep, scores_2_3_p)
          
          # ==== EPITHELIAL SCORE
          scores_0_1_e = full_data_comparison_scores_0_1$epithelial_score[time_ss_0_1_e:t_max_ind]
          scores_0_2_e = full_data_comparison_scores_0_2$epithelial_score[time_ss_0_2_e:t_max_ind]
          scores_0_3_e = full_data_comparison_scores_0_3$epithelial_score[time_ss_0_3_e:t_max_ind]
          scores_1_1_e = full_data_comparison_scores_1_1$epithelial_score[time_ss_1_1_e:t_max_ind]
          scores_1_2_e = full_data_comparison_scores_1_2$epithelial_score[time_ss_1_2_e:t_max_ind]
          scores_1_3_e = full_data_comparison_scores_1_3$epithelial_score[time_ss_1_3_e:t_max_ind]
          scores_2_1_e = full_data_comparison_scores_2_1$epithelial_score[time_ss_2_1_e:t_max_ind]
          scores_2_2_e = full_data_comparison_scores_2_2$epithelial_score[time_ss_2_2_e:t_max_ind]
          scores_2_3_e = full_data_comparison_scores_2_3$epithelial_score[time_ss_2_3_e:t_max_ind]
          
          # Accumulate scores
          scores_0_1_e_keep = c(scores_0_1_e_keep, scores_0_1_e)
          scores_0_2_e_keep = c(scores_0_2_e_keep, scores_0_2_e)
          scores_0_3_e_keep = c(scores_0_3_e_keep, scores_0_3_e)
          scores_1_1_e_keep = c(scores_1_1_e_keep, scores_1_1_e)
          scores_1_2_e_keep = c(scores_1_2_e_keep, scores_1_2_e)
          scores_1_3_e_keep = c(scores_1_3_e_keep, scores_1_3_e)
          scores_2_1_e_keep = c(scores_2_1_e_keep, scores_2_1_e)
          scores_2_2_e_keep = c(scores_2_2_e_keep, scores_2_2_e)
          scores_2_3_e_keep = c(scores_2_3_e_keep, scores_2_3_e)
          
          # Compute oscillation metrics for each signal

          osc_0_1_e = compute_oscillation_metrics(scores_0_1_e)
          osc_0_2_e = compute_oscillation_metrics(scores_0_2_e)
          osc_0_3_e = compute_oscillation_metrics(scores_0_3_e)
          osc_1_1_e = compute_oscillation_metrics(scores_1_1_e)
          osc_1_2_e = compute_oscillation_metrics(scores_1_2_e)
          osc_1_3_e = compute_oscillation_metrics(scores_1_3_e)
          osc_2_1_e = compute_oscillation_metrics(scores_2_1_e)
          osc_2_2_e = compute_oscillation_metrics(scores_2_2_e)
          osc_2_3_e = compute_oscillation_metrics(scores_2_3_e)
          
          osc_0_1_p = compute_oscillation_metrics(scores_0_1_p)
          osc_0_2_p = compute_oscillation_metrics(scores_0_2_p)
          osc_0_3_p = compute_oscillation_metrics(scores_0_3_p)
          osc_1_1_p = compute_oscillation_metrics(scores_1_1_p)
          osc_1_2_p = compute_oscillation_metrics(scores_1_2_p)
          osc_1_3_p = compute_oscillation_metrics(scores_1_3_p)
          osc_2_1_p = compute_oscillation_metrics(scores_2_1_p)
          osc_2_2_p = compute_oscillation_metrics(scores_2_2_p)
          osc_2_3_p = compute_oscillation_metrics(scores_2_3_p)
          
          # --- Tabulate all comparisons ---
          comparison_results = data.frame(
            param_set_id = i,
            replicate_id = rep,
            injury_type  = rep("pathogenic", 18),
            score_type  = c(rep("epithelium", 9), rep("pathogen", 9)),
            macspec_on   = rep(0, 18),
            tregs_on     = rep(0, 18),
            tregs_rnd    = rep(0, 18),
            ros_level    = c(0, 0, 0, 1, 1, 1, 2, 2, 2,
                             0, 0, 0, 1, 1, 1, 2, 2, 2),
            pat_level    = c(1, 2, 3, 1, 2, 3, 1, 2, 3,
                             1, 2, 3, 1, 2, 3, 1, 2, 3),
            ss_start     = c(time_ss_0_1_e, time_ss_0_2_e, time_ss_0_3_e,
                             time_ss_1_1_e, time_ss_1_2_e, time_ss_1_3_e,
                             time_ss_2_1_e, time_ss_2_2_e, time_ss_2_3_e,
                             time_ss_0_1_p, time_ss_0_2_p, time_ss_0_3_p,
                             time_ss_1_1_p, time_ss_1_2_p, time_ss_1_3_p,
                             time_ss_2_1_p, time_ss_2_2_p, time_ss_2_3_p),
            mean_score   = c(mean(scores_0_1_e), mean(scores_0_2_e), mean(scores_0_3_e),
                             mean(scores_1_1_e), mean(scores_1_2_e), mean(scores_1_3_e),
                             mean(scores_2_1_e), mean(scores_2_2_e), mean(scores_2_3_e),
                             mean(scores_0_1_p), mean(scores_0_2_p), mean(scores_0_3_p),
                             mean(scores_1_1_p), mean(scores_1_2_p), mean(scores_1_3_p),
                             mean(scores_2_1_p), mean(scores_2_2_p), mean(scores_2_3_p)),
            sd_score     = c(sd(scores_0_1_e), sd(scores_0_2_e), sd(scores_0_3_e),
                             sd(scores_1_1_e), sd(scores_1_2_e), sd(scores_1_3_e),
                             sd(scores_2_1_e), sd(scores_2_2_e), sd(scores_2_3_e),
                             sd(scores_0_1_p), sd(scores_0_2_p), sd(scores_0_3_p),
                             sd(scores_1_1_p), sd(scores_1_2_p), sd(scores_1_3_p),
                             sd(scores_2_1_p), sd(scores_2_2_p), sd(scores_2_3_p))
            # Oscillation metrics
          )
          
          # # ==== interpretation
          # # ==== acf_peak_val > 0.3–0.4 typically indicates clear oscillations (like your 39518 example)
          # # ==== acf_peak_lag gives you the approximate period
          # # ==== spec_conc > 0.1–0.15 suggests concentrated periodic power vs diffuse noise
          
          # Append to global results
          all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
        }
        else{
          processed_indices = c(processed_indices, i) #add and skip
          message("Skipped one, getting out of loop")
          break
        }
      }
      
      if(dim(all_comparison_results_reps)[1]>0){
        
        # JS Divergence - ALL POSSIBLE COMBINATIONS
        # All comparisons between all 9 conditions with explicit naming
        ros_vals = c(0, 1, 2)
        pat_vals = c(1, 2, 3)
        
        # ==== EPITHELIAL SCORE
        for (r1 in ros_vals) {
          for (p1 in pat_vals) {
            for (r2 in ros_vals) {
              for (p2 in pat_vals) {
                if (r1 != r2 || p1 != p2) {
                  col_name = paste0("d_ros", r1, "_pat", p1, "_treg0_vs_ros", r2, "_pat", p2, "_treg0_e")
                  scores1 = get(paste0("scores_", r1, "_", p1, "_e_keep"))
                  scores2 = get(paste0("scores_", r2, "_", p2, "_e_keep"))
                  all_comparison_results_reps[[col_name]] = calculate_js_divergence(scores1, scores2, method = "fd")[1]
                }
              }
            }
          }
        }
        
        # ==== PATHOGEN ABUNDANCE
        for (r1 in ros_vals) {
          for (p1 in pat_vals) {
            for (r2 in ros_vals) {
              for (p2 in pat_vals) {
                if (r1 != r2 || p1 != p2) {
                  col_name = paste0("d_ros", r1, "_pat", p1, "_treg0_vs_ros", r2, "_pat", p2, "_treg0_p")
                  scores1 = get(paste0("scores_", r1, "_", p1, "_p_keep"))
                  scores2 = get(paste0("scores_", r2, "_", p2, "_p_keep"))
                  all_comparison_results_reps[[col_name]] = calculate_js_divergence(scores1, scores2, method = "fd")[1]
                }
              }
            }
          }
        }
        
        # ====== Mean scores - epithelial score
        all_comparison_results_reps$mean_ros0_pat1_treg0_e = mean(scores_0_1_e_keep)
        all_comparison_results_reps$mean_ros1_pat1_treg0_e = mean(scores_1_1_e_keep)
        all_comparison_results_reps$mean_ros2_pat1_treg0_e = mean(scores_2_1_e_keep)
        all_comparison_results_reps$mean_ros0_pat2_treg0_e = mean(scores_0_2_e_keep)
        all_comparison_results_reps$mean_ros1_pat2_treg0_e = mean(scores_1_2_e_keep)
        all_comparison_results_reps$mean_ros2_pat2_treg0_e = mean(scores_2_2_e_keep)
        all_comparison_results_reps$mean_ros0_pat3_treg0_e = mean(scores_0_3_e_keep)
        all_comparison_results_reps$mean_ros1_pat3_treg0_e = mean(scores_1_3_e_keep)
        all_comparison_results_reps$mean_ros2_pat3_treg0_e = mean(scores_2_3_e_keep)
        
        
        # ====== Mean scores - pathogen abundance
        all_comparison_results_reps$mean_ros0_pat1_treg0_p = mean(scores_0_1_p_keep)
        all_comparison_results_reps$mean_ros1_pat1_treg0_p = mean(scores_1_1_p_keep)
        all_comparison_results_reps$mean_ros2_pat1_treg0_p = mean(scores_2_1_p_keep)
        all_comparison_results_reps$mean_ros0_pat2_treg0_p = mean(scores_0_2_p_keep)
        all_comparison_results_reps$mean_ros1_pat2_treg0_p = mean(scores_1_2_p_keep)
        all_comparison_results_reps$mean_ros2_pat2_treg0_p = mean(scores_2_2_p_keep)
        all_comparison_results_reps$mean_ros0_pat3_treg0_p = mean(scores_0_3_p_keep)
        all_comparison_results_reps$mean_ros1_pat3_treg0_p = mean(scores_1_3_p_keep)
        all_comparison_results_reps$mean_ros2_pat3_treg0_p = mean(scores_2_3_p_keep)
        
        all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
        processed_indices      = c(processed_indices, i)
      }
    }
    # Save after every 10 parameter sets (if total is > 10) or save all (if total is <= 10)
    # if((length(inds2read) > 10 && i_idx %% 10 == 0) || i_idx==length(loop_over)){
    
    if((length(inds2read) > (3*n2) && i_idx %% (3*n2) == 0) || i_idx==length(loop_over)){ #3*n2 to break sync between multiple cores
      
      message("Saving intermediate results after ", i_idx, " parameter sets...")
      
      # Update the list of read indices
      if(!file.exists('./ids_read_abm_patros.rds')){
        saveRDS(processed_indices, './ids_read_abm_patros.rds')
      }else{
        inds_read_old = readRDS('./ids_read_abm_patros.rds')
        inds_read_updated = c(inds_read_old, processed_indices)
        saveRDS(inds_read_updated, './ids_read_abm_patros.rds')
      }
      
      # Update the results data
      if(!file.exists('./data_cpp_read_abm_patros.rds')){
        saveRDS(all_comparison_results, './data_cpp_read_abm_patros.rds')
      }else{
        all_comparison_results_old = readRDS('./data_cpp_read_abm_patros.rds')
        all_comparison_results_combined = rbind(all_comparison_results_old, all_comparison_results)
        saveRDS(all_comparison_results_combined, './data_cpp_read_abm_patros.rds')
      }
      
      # Reset for next batch
      all_comparison_results = data.frame()
      processed_indices = c()
      
      message(paste0("Intermediate save complete with ",3*n2," more param_ids."))
    }
  }
  
  # Final save for any remaining results
  message("Saving final results. Total param_sets processed: ", length(inds2read))
  
  if(!file.exists('./ids_read_abm_patros.rds')){
    saveRDS(processed_indices, './ids_read_abm_patros.rds')
  }else{
    inds_read_old = readRDS('./ids_read_abm_patros.rds')
    inds_read_final = c(inds_read_old, processed_indices)
    saveRDS(inds_read_final, './ids_read_abm_patros.rds')
  }
  
  if(!file.exists('./data_cpp_read_abm_patros.rds')){
    saveRDS(all_comparison_results, './data_cpp_read_abm_patros.rds')
  }else{
    all_comparison_results_old = readRDS('./data_cpp_read_abm_patros.rds')
    all_comparison_results_final = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results_final, './data_cpp_read_abm_patros.rds')
  }
  
  message("Saved.")
}else{
  message("No new pts added.")
}



