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
# Control
files_1_0_0_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds$", full.names = TRUE)
files_1_1_0_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_1_sterile_1_macspec_0_tregs_0_trnd_0.rds$", full.names = TRUE)
# === STERILE
# Test (with ROS)
files_0_0_0_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds$", full.names = TRUE)
files_0_0_0_1_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_0_tregs_1_trnd_0.rds$", full.names = TRUE)
files_0_0_0_1_1 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_0_tregs_1_trnd_1.rds$", full.names = TRUE)
# Test (with ROS), perfect macrophage
files_0_0_1_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_1_tregs_0_trnd_0.rds$", full.names = TRUE)
# === PATHOGENIC
# Test (with ROS)
files_0_1_0_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_1_macspec_0_tregs_0_trnd_0.rds$", full.names = TRUE)
files_0_1_0_1_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_1_macspec_0_tregs_1_trnd_0.rds$", full.names = TRUE)
files_0_1_0_1_1 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_1_macspec_0_tregs_1_trnd_1.rds$", full.names = TRUE)
# Test (with ROS), perfect macrophage
files_0_1_1_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_1_macspec_1_tregs_0_trnd_0.rds$", full.names = TRUE)

indices_1_0_0_0_0 = str_extract(basename(files_1_0_0_0_0), "\\d+") |> as.numeric()
indices_1_1_0_0_0 = str_extract(basename(files_1_1_0_0_0), "\\d+") |> as.numeric()

indices_0_0_0_0_0 = str_extract(basename(files_0_0_0_0_0), "\\d+") |> as.numeric()
indices_0_0_0_1_0 = str_extract(basename(files_0_0_0_1_0), "\\d+") |> as.numeric()
indices_0_0_0_1_1 = str_extract(basename(files_0_0_0_1_1), "\\d+") |> as.numeric()
indices_0_0_1_0_0 = str_extract(basename(files_0_0_1_0_0), "\\d+") |> as.numeric()

indices_0_1_0_0_0 = str_extract(basename(files_0_1_0_0_0), "\\d+") |> as.numeric()
indices_0_1_0_1_0 = str_extract(basename(files_0_1_0_1_0), "\\d+") |> as.numeric()
indices_0_1_0_1_1 = str_extract(basename(files_0_1_0_1_1), "\\d+") |> as.numeric()
indices_0_1_1_0_0 = str_extract(basename(files_0_1_1_0_0), "\\d+") |> as.numeric()

indices = Reduce(intersect, list(
  indices_1_0_0_0_0, indices_1_1_0_0_0,
  indices_0_0_0_0_0, indices_0_0_0_1_0,
  indices_0_0_0_1_1, indices_0_0_1_0_0,
  indices_0_1_0_0_0, indices_0_1_0_1_0,
  indices_0_1_0_1_1, indices_0_1_1_0_0
))

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('./ids_read_abm.rds')){
  inds_read = readRDS('./ids_read_abm.rds')
}else{
  inds_read = c()
}

inds2read = sort(setdiff(indices,inds_read))
length(inds2read)

if(length(inds2read)>0){
  # Track indices that have been successfully processed
  processed_indices = c()
  
  for (i_idx in seq_along(inds2read)){
    
    i = inds2read[i_idx]
    message("Processing param_set_", i)
    # Check file sizes for this parameter set
    files_to_check = c(
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_1_trnd_1.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_1_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_1_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_1_trnd_1.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_1_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_1_macspec_0_tregs_0_trnd_0.rds')
    )
    if(any(file.info(files_to_check)$size<100000)){
      processed_indices      = c(processed_indices, i) #add and skip
      message("Skipped one")
    }else{
      # Control files
      results_1_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
      results_1_1_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_1_macspec_0_tregs_0_trnd_0.rds'))
      
      # Pathogenic (sterile_0) test files
      results_0_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
      results_0_0_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_1_trnd_0.rds'))
      results_0_0_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_1_trnd_1.rds'))
      results_0_0_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_1_tregs_0_trnd_0.rds'))
      
      # Sterile (sterile_1) test files
      results_0_1_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_0_trnd_0.rds'))
      results_0_1_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_1_trnd_0.rds'))
      results_0_1_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_0_tregs_1_trnd_1.rds'))
      results_0_1_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_1_macspec_1_tregs_0_trnd_0.rds'))
      
      results = rbind(
        results_1_0_0_0_0, results_1_1_0_0_0,
        results_0_0_0_0_0, results_0_0_0_1_0,
        results_0_0_0_1_1, results_0_0_1_0_0,
        results_0_1_0_0_0, results_0_1_0_1_0,
        results_0_1_0_1_1, results_0_1_1_0_0
      )
      
      full_data_comparison = results %>% dplyr::select(param_set_id, control, sterile, macspec_on, tregs_on, randomize_tregs, rep_id, t, time_ss, epithelial_score, pathogen)
      min_reps  = min(full_data_comparison$rep_id)
      max_reps  = max(full_data_comparison$rep_id)
      t_max_ind = max(full_data_comparison$t)
      
      # ====== PATHOGEN
      # Control
      scores_1_0_0_0_0_p_keep = c()  # control_1_sterile_0_macspec_0_tregs_0_trnd_0
      scores_1_1_0_0_0_p_keep = c()  # control_1_sterile_1_macspec_0_tregs_0_trnd_0
      
      # Pathogenic (sterile_0) test files
      scores_0_0_0_0_0_p_keep = c()  # control_0_sterile_0_macspec_0_tregs_0_trnd_0
      scores_0_0_0_1_0_p_keep = c()  # control_0_sterile_0_macspec_0_tregs_1_trnd_0
      scores_0_0_0_1_1_p_keep = c()  # control_0_sterile_0_macspec_0_tregs_1_trnd_1
      scores_0_0_1_0_0_p_keep = c()  # control_0_sterile_0_macspec_1_tregs_0_trnd_0
      
      # Sterile (sterile_1) test files
      scores_0_1_0_0_0_p_keep = c()  # control_0_sterile_1_macspec_0_tregs_0_trnd_0
      scores_0_1_0_1_0_p_keep = c()  # control_0_sterile_1_macspec_0_tregs_1_trnd_0
      scores_0_1_0_1_1_p_keep = c()  # control_0_sterile_1_macspec_0_tregs_1_trnd_1
      scores_0_1_1_0_0_p_keep = c()  # control_0_sterile_1_macspec_1_tregs_0_trnd_0
      
      # ===== EPITHELIAL SCORE
      # Control
      scores_1_0_0_0_0_e_keep = c()  # control_1_sterile_0_macspec_0_tregs_0_trnd_0
      scores_1_1_0_0_0_e_keep = c()  # control_1_sterile_1_macspec_0_tregs_0_trnd_0
      
      # Pathogenic (sterile_0) test files
      scores_0_0_0_0_0_e_keep = c()  # control_0_sterile_0_macspec_0_tregs_0_trnd_0
      scores_0_0_0_1_0_e_keep = c()  # control_0_sterile_0_macspec_0_tregs_1_trnd_0
      scores_0_0_0_1_1_e_keep = c()  # control_0_sterile_0_macspec_0_tregs_1_trnd_1
      scores_0_0_1_0_0_e_keep = c()  # control_0_sterile_0_macspec_1_tregs_0_trnd_0
      
      # Sterile (sterile_1) test files
      scores_0_1_0_0_0_e_keep = c()  # control_0_sterile_1_macspec_0_tregs_0_trnd_0
      scores_0_1_0_1_0_e_keep = c()  # control_0_sterile_1_macspec_0_tregs_1_trnd_0
      scores_0_1_0_1_1_e_keep = c()  # control_0_sterile_1_macspec_0_tregs_1_trnd_1
      scores_0_1_1_0_0_e_keep = c()  # control_0_sterile_1_macspec_1_tregs_0_trnd_0
      
      all_comparison_results_reps = data.frame()
      
      for (rep in min_reps:max_reps) {
        
        #### CONTROL
        # sterile_0 (no pathogen)
        full_data_comparison_scores_1_0_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==1 & sterile==0 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        # sterile_1 (with pathogen)
        full_data_comparison_scores_1_1_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==1 & sterile==1 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        
        #### Pathogenic (sterile_0) - Test with ROS
        # tregs OFF
        full_data_comparison_scores_0_0_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        # tregs ON
        full_data_comparison_scores_0_0_0_1_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==0 & tregs_on==1 & randomize_tregs==0)
        # tregs ON, BUT ARE random
        full_data_comparison_scores_0_0_0_1_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==0 & tregs_on==1 & randomize_tregs==1)
        # perfect macrophage, tregs OFF
        full_data_comparison_scores_0_0_1_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==1 & tregs_on==0 & randomize_tregs==0)
        
        #### Sterile (sterile_1) - Test with ROS
        # tregs OFF
        full_data_comparison_scores_0_1_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==1 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        # tregs ON
        full_data_comparison_scores_0_1_0_1_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==1 & macspec_on==0 & tregs_on==1 & randomize_tregs==0)
        # tregs ON, BUT ARE random
        full_data_comparison_scores_0_1_0_1_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==1 & macspec_on==0 & tregs_on==1 & randomize_tregs==1)
        # perfect macrophage, tregs OFF
        full_data_comparison_scores_0_1_1_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==1 & macspec_on==1 & tregs_on==0 & randomize_tregs==0)
        
        # --- Steady-state detection ---
        time_ss_1_0_0_0_0 = unique(full_data_comparison_scores_1_0_0_0_0$time_ss)
        time_ss_1_1_0_0_0 = unique(full_data_comparison_scores_1_1_0_0_0$time_ss)
        
        time_ss_0_0_0_0_0 = unique(full_data_comparison_scores_0_0_0_0_0$time_ss)
        time_ss_0_0_0_1_0 = unique(full_data_comparison_scores_0_0_0_1_0$time_ss)
        time_ss_0_0_0_1_1 = unique(full_data_comparison_scores_0_0_0_1_1$time_ss)
        time_ss_0_0_1_0_0 = unique(full_data_comparison_scores_0_0_1_0_0$time_ss)
        
        time_ss_0_1_0_0_0 = unique(full_data_comparison_scores_0_1_0_0_0$time_ss)
        time_ss_0_1_0_1_0 = unique(full_data_comparison_scores_0_1_0_1_0$time_ss)
        time_ss_0_1_0_1_1 = unique(full_data_comparison_scores_0_1_0_1_1$time_ss)
        time_ss_0_1_1_0_0 = unique(full_data_comparison_scores_0_1_1_0_0$time_ss)
        
        time_ss_vec = c(
          time_ss_1_0_0_0_0, time_ss_1_1_0_0_0,
          time_ss_0_0_0_0_0, time_ss_0_0_0_1_0, 
          time_ss_0_0_0_1_1, time_ss_0_0_1_0_0,
          time_ss_0_1_0_0_0, time_ss_0_1_0_1_0, 
          time_ss_0_1_0_1_1, time_ss_0_1_1_0_0
        )
        
        if(!any(is.na(time_ss_vec))){
          
          # ==== PATHOGEN ABUNDANCE
          # Control
          scores_1_0_0_0_0_p = full_data_comparison_scores_1_0_0_0_0$pathogen[time_ss_1_0_0_0_0:t_max_ind]
          scores_1_1_0_0_0_p = full_data_comparison_scores_1_1_0_0_0$pathogen[time_ss_1_1_0_0_0:t_max_ind]
          
          # Pathogenic (sterile_0)
          scores_0_0_0_0_0_p = full_data_comparison_scores_0_0_0_0_0$pathogen[time_ss_0_0_0_0_0:t_max_ind]
          scores_0_0_0_1_0_p = full_data_comparison_scores_0_0_0_1_0$pathogen[time_ss_0_0_0_1_0:t_max_ind]
          scores_0_0_0_1_1_p = full_data_comparison_scores_0_0_0_1_1$pathogen[time_ss_0_0_0_1_1:t_max_ind]
          scores_0_0_1_0_0_p = full_data_comparison_scores_0_0_1_0_0$pathogen[time_ss_0_0_1_0_0:t_max_ind]
          
          # Sterile (sterile_1)
          scores_0_1_0_0_0_p = full_data_comparison_scores_0_1_0_0_0$pathogen[time_ss_0_1_0_0_0:t_max_ind]
          scores_0_1_0_1_0_p = full_data_comparison_scores_0_1_0_1_0$pathogen[time_ss_0_1_0_1_0:t_max_ind]
          scores_0_1_0_1_1_p = full_data_comparison_scores_0_1_0_1_1$pathogen[time_ss_0_1_0_1_1:t_max_ind]
          scores_0_1_1_0_0_p = full_data_comparison_scores_0_1_1_0_0$pathogen[time_ss_0_1_1_0_0:t_max_ind]
          
          # Accumulate scores
          scores_1_0_0_0_0_p_keep = c(scores_1_0_0_0_0_p_keep, scores_1_0_0_0_0_p)
          scores_1_1_0_0_0_p_keep = c(scores_1_1_0_0_0_p_keep, scores_1_1_0_0_0_p)
          
          scores_0_0_0_0_0_p_keep = c(scores_0_0_0_0_0_p_keep, scores_0_0_0_0_0_p)
          scores_0_0_0_1_0_p_keep = c(scores_0_0_0_1_0_p_keep, scores_0_0_0_1_0_p)
          scores_0_0_0_1_1_p_keep = c(scores_0_0_0_1_1_p_keep, scores_0_0_0_1_1_p)
          scores_0_0_1_0_0_p_keep = c(scores_0_0_1_0_0_p_keep, scores_0_0_1_0_0_p)
          
          scores_0_1_0_0_0_p_keep = c(scores_0_1_0_0_0_p_keep, scores_0_1_0_0_0_p)
          scores_0_1_0_1_0_p_keep = c(scores_0_1_0_1_0_p_keep, scores_0_1_0_1_0_p)
          scores_0_1_0_1_1_p_keep = c(scores_0_1_0_1_1_p_keep, scores_0_1_0_1_1_p)
          scores_0_1_1_0_0_p_keep = c(scores_0_1_1_0_0_p_keep, scores_0_1_1_0_0_p)
          
          
          # ==== EPITHELIAL SCORE
          # Control
          scores_1_0_0_0_0_e = full_data_comparison_scores_1_0_0_0_0$epithelial_score[time_ss_1_0_0_0_0:t_max_ind]
          scores_1_1_0_0_0_e = full_data_comparison_scores_1_1_0_0_0$epithelial_score[time_ss_1_1_0_0_0:t_max_ind]
          
          # Pathogenic (sterile_0)
          scores_0_0_0_0_0_e = full_data_comparison_scores_0_0_0_0_0$epithelial_score[time_ss_0_0_0_0_0:t_max_ind]
          scores_0_0_0_1_0_e = full_data_comparison_scores_0_0_0_1_0$epithelial_score[time_ss_0_0_0_1_0:t_max_ind]
          scores_0_0_0_1_1_e = full_data_comparison_scores_0_0_0_1_1$epithelial_score[time_ss_0_0_0_1_1:t_max_ind]
          scores_0_0_1_0_0_e = full_data_comparison_scores_0_0_1_0_0$epithelial_score[time_ss_0_0_1_0_0:t_max_ind]
          
          # Sterile (sterile_1)
          scores_0_1_0_0_0_e = full_data_comparison_scores_0_1_0_0_0$epithelial_score[time_ss_0_1_0_0_0:t_max_ind]
          scores_0_1_0_1_0_e = full_data_comparison_scores_0_1_0_1_0$epithelial_score[time_ss_0_1_0_1_0:t_max_ind]
          scores_0_1_0_1_1_e = full_data_comparison_scores_0_1_0_1_1$epithelial_score[time_ss_0_1_0_1_1:t_max_ind]
          scores_0_1_1_0_0_e = full_data_comparison_scores_0_1_1_0_0$epithelial_score[time_ss_0_1_1_0_0:t_max_ind]
          
          # Accumulate scores
          scores_1_0_0_0_0_e_keep = c(scores_1_0_0_0_0_e_keep, scores_1_0_0_0_0_e)
          scores_1_1_0_0_0_e_keep = c(scores_1_1_0_0_0_e_keep, scores_1_1_0_0_0_e)
          
          scores_0_0_0_0_0_e_keep = c(scores_0_0_0_0_0_e_keep, scores_0_0_0_0_0_e)
          scores_0_0_0_1_0_e_keep = c(scores_0_0_0_1_0_e_keep, scores_0_0_0_1_0_e)
          scores_0_0_0_1_1_e_keep = c(scores_0_0_0_1_1_e_keep, scores_0_0_0_1_1_e)
          scores_0_0_1_0_0_e_keep = c(scores_0_0_1_0_0_e_keep, scores_0_0_1_0_0_e)
          
          scores_0_1_0_0_0_e_keep = c(scores_0_1_0_0_0_e_keep, scores_0_1_0_0_0_e)
          scores_0_1_0_1_0_e_keep = c(scores_0_1_0_1_0_e_keep, scores_0_1_0_1_0_e)
          scores_0_1_0_1_1_e_keep = c(scores_0_1_0_1_1_e_keep, scores_0_1_0_1_1_e)
          scores_0_1_1_0_0_e_keep = c(scores_0_1_1_0_0_e_keep, scores_0_1_1_0_0_e)
          
          # --- Tabulate all comparisons ---
          comparison_results = data.frame(
            param_set_id = i,
            replicate_id = rep,
            control      = c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
            injury_type  = c("pathogenic", "sterile", 
                             "pathogenic", "pathogenic", 
                             "pathogenic", "pathogenic", 
                             "sterile", "sterile", 
                             "sterile", "sterile",
                             "pathogenic", "sterile", 
                             "pathogenic", "pathogenic", 
                             "pathogenic", "pathogenic", 
                             "sterile", "sterile", 
                             "sterile", "sterile"),
            macspec_on   = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
            tregs_on     = c(0, 0, 0, 1, 1, 0, 0, 1, 1, 0,
                             0, 0, 0, 1, 1, 0, 0, 1, 1, 0),
            tregs_rnd    = c(0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 1, 0),
            ss_start     = c(time_ss_1_0_0_0_0, time_ss_1_1_0_0_0, 
                             time_ss_0_0_0_0_0, time_ss_0_0_0_1_0, 
                             time_ss_0_0_0_1_1, time_ss_0_0_1_0_0,
                             time_ss_0_1_0_0_0, time_ss_0_1_0_1_0, 
                             time_ss_0_1_0_1_1, time_ss_0_1_1_0_0,
                             time_ss_1_0_0_0_0, time_ss_1_1_0_0_0, 
                             time_ss_0_0_0_0_0, time_ss_0_0_0_1_0, 
                             time_ss_0_0_0_1_1, time_ss_0_0_1_0_0,
                             time_ss_0_1_0_0_0, time_ss_0_1_0_1_0, 
                             time_ss_0_1_0_1_1, time_ss_0_1_1_0_0),
            mean_score   = c(mean(scores_1_0_0_0_0_e), mean(scores_1_1_0_0_0_e),
                             mean(scores_0_0_0_0_0_e), mean(scores_0_0_0_1_0_e), 
                             mean(scores_0_0_0_1_1_e), mean(scores_0_0_1_0_0_e),
                             mean(scores_0_1_0_0_0_e), mean(scores_0_1_0_1_0_e), 
                             mean(scores_0_1_0_1_1_e), mean(scores_0_1_1_0_0_e),
                             mean(scores_1_0_0_0_0_p), mean(scores_1_1_0_0_0_p),
                             mean(scores_0_0_0_0_0_p), mean(scores_0_0_0_1_0_p), 
                             mean(scores_0_0_0_1_1_p), mean(scores_0_0_1_0_0_p),
                             mean(scores_0_1_0_0_0_p), mean(scores_0_1_0_1_0_p), 
                             mean(scores_0_1_0_1_1_p), mean(scores_0_1_1_0_0_p))
          )
          
          # Append to global results
          all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
        }
      }
      
      # ====== EPITHELIAL SCORE
      # JS Divergence comparisons
      # Control vs tregs OFF (baseline effect of ROS/injury)
      all_comparison_results_reps$d_ctrl_vs_test_sterile_e   = calculate_js_divergence(scores_1_1_0_0_0_e_keep, scores_0_1_0_0_0_e_keep, method = "fd")[1]
      all_comparison_results_reps$d_ctrl_vs_test_pathogen_e  = calculate_js_divergence(scores_1_0_0_0_0_e_keep, scores_0_0_0_0_0_e_keep, method = "fd")[1]
      
      # Tregs ON vs OFF (effect of tregs)
      all_comparison_results_reps$d_tregs_on_vs_off_sterile_e  = calculate_js_divergence(scores_0_1_0_1_0_e_keep, scores_0_1_0_0_0_e_keep, method = "fd")[1]
      all_comparison_results_reps$d_tregs_on_vs_off_pathogen_e = calculate_js_divergence(scores_0_0_0_1_0_e_keep, scores_0_0_0_0_0_e_keep, method = "fd")[1]
      
      # Tregs random vs targeted (specificity of tregs)
      all_comparison_results_reps$d_tregs_rnd_vs_nrnd_sterile_e  = calculate_js_divergence(scores_0_1_0_1_1_e_keep, scores_0_1_0_1_0_e_keep, method = "fd")[1]
      all_comparison_results_reps$d_tregs_rnd_vs_nrnd_pathogen_e = calculate_js_divergence(scores_0_0_0_1_1_e_keep, scores_0_0_0_1_0_e_keep, method = "fd")[1]
      
      # Perfect macrophage vs tregs OFF (macrophage specificity effect)
      all_comparison_results_reps$d_macspec_vs_tregs_off_sterile_e  = calculate_js_divergence(scores_0_1_1_0_0_e_keep, scores_0_1_0_0_0_e_keep, method = "fd")[1]
      all_comparison_results_reps$d_macspec_vs_tregs_off_pathogen_e = calculate_js_divergence(scores_0_0_1_0_0_e_keep, scores_0_0_0_0_0_e_keep, method = "fd")[1]
      
      # Perfect macrophage vs tregs ON (comparing two rescue mechanisms)
      all_comparison_results_reps$d_macspec_vs_tregs_on_sterile_e  = calculate_js_divergence(scores_0_1_1_0_0_e_keep, scores_0_1_0_1_0_e_keep, method = "fd")[1]
      all_comparison_results_reps$d_macspec_vs_tregs_on_pathogen_e = calculate_js_divergence(scores_0_0_1_0_0_e_keep, scores_0_0_0_1_0_e_keep, method = "fd")[1]
      
      # ===== PATHOGEN
      # JS Divergence comparisons
      # Control vs tregs OFF (baseline effect of ROS/injury)
      all_comparison_results_reps$d_ctrl_vs_test_sterile_p   = calculate_js_divergence(scores_1_1_0_0_0_p_keep, scores_0_1_0_0_0_p_keep, method = "fd")[1]
      all_comparison_results_reps$d_ctrl_vs_test_pathogen_p  = calculate_js_divergence(scores_1_0_0_0_0_p_keep, scores_0_0_0_0_0_p_keep, method = "fd")[1]
      
      # Tregs ON vs OFF (effect of tregs)
      all_comparison_results_reps$d_tregs_on_vs_off_sterile_p  = calculate_js_divergence(scores_0_1_0_1_0_p_keep, scores_0_1_0_0_0_p_keep, method = "fd")[1]
      all_comparison_results_reps$d_tregs_on_vs_off_pathogen_p = calculate_js_divergence(scores_0_0_0_1_0_p_keep, scores_0_0_0_0_0_p_keep, method = "fd")[1]
      
      # Tregs random vs targeted (specificity of tregs)
      all_comparison_results_reps$d_tregs_rnd_vs_nrnd_sterile_p  = calculate_js_divergence(scores_0_1_0_1_1_p_keep, scores_0_1_0_1_0_p_keep, method = "fd")[1]
      all_comparison_results_reps$d_tregs_rnd_vs_nrnd_pathogen_p = calculate_js_divergence(scores_0_0_0_1_1_p_keep, scores_0_0_0_1_0_p_keep, method = "fd")[1]
      
      # Perfect macrophage vs tregs OFF (macrophage specificity effect)
      all_comparison_results_reps$d_macspec_vs_tregs_off_sterile_p       = calculate_js_divergence(scores_0_1_1_0_0_p_keep, scores_0_1_0_0_0_p_keep, method = "fd")[1]
      all_comparison_results_reps$d_macspec_vs_tregs_off_pathogen_p = calculate_js_divergence(scores_0_0_1_0_0_p_keep, scores_0_0_0_0_0_p_keep, method = "fd")[1]
      
      # Perfect macrophage vs tregs ON (comparing two rescue mechanisms)
      all_comparison_results_reps$d_macspec_vs_tregs_on_sterile_p  = calculate_js_divergence(scores_0_1_1_0_0_p_keep, scores_0_1_0_1_0_p_keep, method = "fd")[1]
      all_comparison_results_reps$d_macspec_vs_tregs_on_pathogen_p = calculate_js_divergence(scores_0_0_1_0_0_p_keep, scores_0_0_0_1_0_p_keep, method = "fd")[1]
      
      
      # Mean scores
      # Control
      all_comparison_results_reps$mean_ctrl_sterile_e  = mean(scores_1_1_0_0_0_e_keep)
      all_comparison_results_reps$mean_ctrl_pathogen_e = mean(scores_1_0_0_0_0_e_keep)
      
      # Test - tregs OFF
      all_comparison_results_reps$mean_tregs_off_sterile_e  = mean(scores_0_1_0_0_0_e_keep)
      all_comparison_results_reps$mean_tregs_off_pathogen_e = mean(scores_0_0_0_0_0_e_keep)
      
      # Test - tregs ON (targeted)
      all_comparison_results_reps$mean_tregs_on_sterile_e  = mean(scores_0_1_0_1_0_e_keep)
      all_comparison_results_reps$mean_tregs_on_pathogen_e = mean(scores_0_0_0_1_0_e_keep)
      
      # Test - tregs ON (random)
      all_comparison_results_reps$mean_tregs_rnd_sterile_e  = mean(scores_0_1_0_1_1_e_keep)
      all_comparison_results_reps$mean_tregs_rnd_pathogen_e = mean(scores_0_0_0_1_1_e_keep)
      
      # Test - perfect macrophage
      all_comparison_results_reps$mean_macspec_sterile_e  = mean(scores_0_1_1_0_0_e_keep)
      all_comparison_results_reps$mean_macspec_pathogen_e = mean(scores_0_0_1_0_0_e_keep)
      
      # SD scores
      # Control
      all_comparison_results_reps$sd_ctrl_sterile_e  = sd(scores_1_1_0_0_0_e_keep)
      all_comparison_results_reps$sd_ctrl_pathogen_e = sd(scores_1_0_0_0_0_e_keep)
      
      # Test - tregs OFF
      all_comparison_results_reps$sd_tregs_off_sterile_e  = sd(scores_0_1_0_0_0_e_keep)
      all_comparison_results_reps$sd_tregs_off_pathogen_e = sd(scores_0_0_0_0_0_e_keep)
      
      # Test - tregs ON (targeted)
      all_comparison_results_reps$sd_tregs_on_sterile_e  = sd(scores_0_1_0_1_0_e_keep)
      all_comparison_results_reps$sd_tregs_on_pathogen_e = sd(scores_0_0_0_1_0_e_keep)
      
      # Test - tregs ON (random)
      all_comparison_results_reps$sd_tregs_rnd_sterile_e  = sd(scores_0_1_0_1_1_e_keep)
      all_comparison_results_reps$sd_tregs_rnd_pathogen_e = sd(scores_0_0_0_1_1_e_keep)
      
      # Test - perfect macrophage
      all_comparison_results_reps$sd_macspec_sterile_e  = sd(scores_0_1_1_0_0_e_keep)
      all_comparison_results_reps$sd_macspec_pathogen_e = sd(scores_0_0_1_0_0_e_keep)
      
      # Mean scores
      # Control
      all_comparison_results_reps$mean_ctrl_sterile_p  = mean(scores_1_1_0_0_0_p_keep)
      all_comparison_results_reps$mean_ctrl_pathogen_p = mean(scores_1_0_0_0_0_p_keep)
      
      # Test - tregs OFF
      all_comparison_results_reps$mean_tregs_off_sterile_p  = mean(scores_0_1_0_0_0_p_keep)
      all_comparison_results_reps$mean_tregs_off_pathogen_p = mean(scores_0_0_0_0_0_p_keep)
      
      # Test - tregs ON (targeted)
      all_comparison_results_reps$mean_tregs_on_sterile_p  = mean(scores_0_1_0_1_0_p_keep)
      all_comparison_results_reps$mean_tregs_on_pathogen_p = mean(scores_0_0_0_1_0_p_keep)
      
      # Test - tregs ON (random)
      all_comparison_results_reps$mean_tregs_rnd_sterile_p  = mean(scores_0_1_0_1_1_p_keep)
      all_comparison_results_reps$mean_tregs_rnd_pathogen_p = mean(scores_0_0_0_1_1_p_keep)
      
      # Test - perfect macrophage
      all_comparison_results_reps$mean_macspec_sterile_p  = mean(scores_0_1_1_0_0_p_keep)
      all_comparison_results_reps$mean_macspec_pathogen_p = mean(scores_0_0_1_0_0_p_keep)
      
      # SD scores
      # Control
      all_comparison_results_reps$sd_ctrl_sterile_p  = sd(scores_1_1_0_0_0_p_keep)
      all_comparison_results_reps$sd_ctrl_pathogen_p = sd(scores_1_0_0_0_0_p_keep)
      
      # Test - tregs OFF
      all_comparison_results_reps$sd_tregs_off_sterile_p  = sd(scores_0_1_0_0_0_p_keep)
      all_comparison_results_reps$sd_tregs_off_pathogen_p = sd(scores_0_0_0_0_0_p_keep)
      
      # Test - tregs ON (targeted)
      all_comparison_results_reps$sd_tregs_on_sterile_p  = sd(scores_0_1_0_1_0_p_keep)
      all_comparison_results_reps$sd_tregs_on_pathogen_p = sd(scores_0_0_0_1_0_p_keep)
      
      # Test - tregs ON (random)
      all_comparison_results_reps$sd_tregs_rnd_sterile_p  = sd(scores_0_1_0_1_1_p_keep)
      all_comparison_results_reps$sd_tregs_rnd_pathogen_p = sd(scores_0_0_0_1_1_p_keep)
      
      # Test - perfect macrophage
      all_comparison_results_reps$sd_macspec_sterile_p  = sd(scores_0_1_1_0_0_p_keep)
      all_comparison_results_reps$sd_macspec_pathogen_p = sd(scores_0_0_1_0_0_p_keep)
      
      all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
      processed_indices      = c(processed_indices, i)
    }
    # Save after every 10 parameter sets (if total is > 10)
    if(length(inds2read) > 10 && i_idx %% 10 == 0){
      message("Saving intermediate results after ", i_idx, " parameter sets...")
      
      # Update the list of read indices
      if(!file.exists('./ids_read_abm.rds')){
        saveRDS(processed_indices, './ids_read_abm.rds')
      }else{
        inds_read_old = readRDS('./ids_read_abm.rds')
        inds_read_updated = c(inds_read_old, processed_indices)
        saveRDS(inds_read_updated, './ids_read_abm.rds')
      }
      
      # Update the results data
      if(!file.exists('./data_cpp_read_abm.rds')){
        saveRDS(all_comparison_results, './data_cpp_read_abm.rds')
      }else{
        all_comparison_results_old = readRDS('./data_cpp_read_abm.rds')
        all_comparison_results_combined = rbind(all_comparison_results_old, all_comparison_results)
        saveRDS(all_comparison_results_combined, './data_cpp_read_abm.rds')
      }
      
      # Reset for next batch
      all_comparison_results = data.frame()
      processed_indices = c()
      
      message("Intermediate save complete with 10 more param_ids.")
    }
  }
  
  # Final save for any remaining results
  message("Saving final results. Total param_sets processed: ", length(inds2read))
  
  if(!file.exists('./ids_read_abm.rds')){
    saveRDS(processed_indices, './ids_read_abm.rds')
  }else{
    inds_read_old = readRDS('./ids_read_abm.rds')
    inds_read_final = c(inds_read_old, processed_indices)
    saveRDS(inds_read_final, './ids_read_abm.rds')
  }
  
  if(!file.exists('./data_cpp_read_abm.rds')){
    saveRDS(all_comparison_results, './data_cpp_read_abm.rds')
  }else{
    all_comparison_results_old = readRDS('./data_cpp_read_abm.rds')
    all_comparison_results_final = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results_final, './data_cpp_read_abm.rds')
  }
  
  # message("Reprinting space.")
  # sd_tol_in   = Inf
  # jsd_th      = 0.3
  # tol_in      = 125*0.2
  # M1_M2_diff  = 1
  # source('./DLL_datazanalyse_abm_regions_with_pathogenic.R')
  # source('./DLL_datazanalyse_abm_all.R')
  
}else{
  message("No new pts added.")
}




