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
# Test (with ROS)
files_0_0_0_0_0 = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds$", full.names = TRUE)

indices_1_0_0_0_0 = str_extract(basename(files_1_0_0_0_0), "\\d+") |> as.numeric()
indices_0_0_0_0_0 = str_extract(basename(files_0_0_0_0_0), "\\d+") |> as.numeric()

indices = Reduce(intersect, list(
  indices_1_0_0_0_0,
  indices_0_0_0_0_0
))

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('./ids_read_abm_ros_vs_ctrl.rds')){
  inds_read = readRDS('./ids_read_abm_ros_vs_ctrl.rds')
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
    # Check file sizes for this parameter set
    files_to_check = c(
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds'),
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds')
    )
    if(any(file.info(files_to_check)$size<100000)){
      processed_indices      = c(processed_indices, i) #add and skip
      message("Skipped one")
    }else{
      # Control files
      results_1_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
      # Pathogenic (sterile_0) test files
      results_0_0_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i, '_control_0_sterile_0_macspec_0_tregs_0_trnd_0.rds'))
      
      results = rbind(
        results_1_0_0_0_0, 
        results_0_0_0_0_0
      )
      
      full_data_comparison = results %>% dplyr::select(param_set_id, control, sterile, macspec_on, tregs_on, randomize_tregs, rep_id, t, time_ss, epithelial_score, pathogen)
      min_reps  = min(full_data_comparison$rep_id)
      max_reps  = max(full_data_comparison$rep_id)
      t_max_ind = max(full_data_comparison$t)
      
      # ====== PATHOGEN
      # Control
      scores_1_0_0_0_0_p_keep = c()  # control_1_sterile_0_macspec_0_tregs_0_trnd_0
      # Pathogenic (sterile_0) test files
      scores_0_0_0_0_0_p_keep = c()  # control_0_sterile_0_macspec_0_tregs_0_trnd_0
      
      # === EPITHELIUM
      scores_1_0_0_0_0_e_keep = c()  # control_1_sterile_0_macspec_0_tregs_0_trnd_0
      # Pathogenic (sterile_0) test files
      scores_0_0_0_0_0_e_keep = c()  # control_0_sterile_0_macspec_0_tregs_0_trnd_0
      
      all_comparison_results_reps = data.frame()
      
      for (rep in min_reps:max_reps) {
        
        #### CONTROL
        # sterile_0 (no pathogen)
        full_data_comparison_scores_1_0_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==1 & sterile==0 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        
        #### Pathogenic (sterile_0) - Test with ROS
        # tregs OFF
        full_data_comparison_scores_0_0_0_0_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==0 & tregs_on==0 & randomize_tregs==0)
        
        # --- Steady-state detection ---
        time_ss_1_0_0_0_0 = unique(full_data_comparison_scores_1_0_0_0_0$time_ss)
        time_ss_0_0_0_0_0 = unique(full_data_comparison_scores_0_0_0_0_0$time_ss)
        
        time_ss_vec = c(
          time_ss_1_0_0_0_0,
          time_ss_0_0_0_0_0
        )
        
        if(!any(is.na(time_ss_vec))){
          
          # ==== PATHOGEN ABUNDANCE
          # Control
          scores_1_0_0_0_0_p = full_data_comparison_scores_1_0_0_0_0$pathogen[time_ss_1_0_0_0_0:t_max_ind]
          # Pathogenic (sterile_0)
          scores_0_0_0_0_0_p = full_data_comparison_scores_0_0_0_0_0$pathogen[time_ss_0_0_0_0_0:t_max_ind]
          # Accumulate scores
          scores_1_0_0_0_0_p_keep = c(scores_1_0_0_0_0_p_keep, scores_1_0_0_0_0_p)
          scores_0_0_0_0_0_p_keep = c(scores_0_0_0_0_0_p_keep, scores_0_0_0_0_0_p)
          
          # ==== EPITHELIAL SCORE
          # Control
          scores_1_0_0_0_0_e = full_data_comparison_scores_1_0_0_0_0$epithelial_score[time_ss_1_0_0_0_0:t_max_ind]
          # Pathogenic (sterile_0)
          scores_0_0_0_0_0_e = full_data_comparison_scores_0_0_0_0_0$epithelial_score[time_ss_0_0_0_0_0:t_max_ind]
          # Accumulate scores
          scores_1_0_0_0_0_e_keep = c(scores_1_0_0_0_0_e_keep, scores_1_0_0_0_0_e)
          scores_0_0_0_0_0_e_keep = c(scores_0_0_0_0_0_e_keep, scores_0_0_0_0_0_e)
          
          # Compute oscillation metrics for each signal
          osc_1_0_0_0_0_e = compute_oscillation_metrics(scores_1_0_0_0_0_e)
          osc_0_0_0_0_0_e = compute_oscillation_metrics(scores_0_0_0_0_0_e)
          osc_1_0_0_0_0_p = compute_oscillation_metrics(scores_1_0_0_0_0_p)
          osc_0_0_0_0_0_p = compute_oscillation_metrics(scores_0_0_0_0_0_p)
          
          # --- Tabulate all comparisons ---
          comparison_results = data.frame(
            param_set_id = i,
            replicate_id = rep,
            control      = c(1, 0,
                             1, 0),
            injury_type  = c("pathogenic", "pathogenic",
                             "pathogenic", "pathogenic"),
            score_type  = c("epithelium", "epithelium",
                            "pathogen", "pathogen"),
            macspec_on   = c(0, 0, 0, 0),
            tregs_on     = c(0, 0, 0, 0),
            tregs_rnd    = c(0, 0, 0, 0),
            ss_start     = c(time_ss_1_0_0_0_0, 
                             time_ss_0_0_0_0_0,
                             time_ss_1_0_0_0_0, 
                             time_ss_0_0_0_0_0),
            mean_score   = c(mean(scores_1_0_0_0_0_e),
                             mean(scores_0_0_0_0_0_e),
                             mean(scores_1_0_0_0_0_p),
                             mean(scores_0_0_0_0_0_p)),
            sd_score   = c(sd(scores_1_0_0_0_0_e),
                           sd(scores_0_0_0_0_0_e),
                           sd(scores_1_0_0_0_0_p),
                           sd(scores_0_0_0_0_0_p)),
            # Oscillation metrics
            acf_peak_lag = c(osc_1_0_0_0_0_e$acf_peak_lag, osc_0_0_0_0_0_e$acf_peak_lag,
                             osc_1_0_0_0_0_p$acf_peak_lag, osc_0_0_0_0_0_p$acf_peak_lag),
            acf_peak_val = c(osc_1_0_0_0_0_e$acf_peak_value, osc_0_0_0_0_0_e$acf_peak_value,
                             osc_1_0_0_0_0_p$acf_peak_value, osc_0_0_0_0_0_p$acf_peak_value),
            spec_conc    = c(osc_1_0_0_0_0_e$spectral_concentration, osc_0_0_0_0_0_e$spectral_concentration,
                             osc_1_0_0_0_0_p$spectral_concentration, osc_0_0_0_0_0_p$spectral_concentration),
            spec_sharp   = c(osc_1_0_0_0_0_e$spectral_sharpness, osc_0_0_0_0_0_e$spectral_sharpness,
                             osc_1_0_0_0_0_p$spectral_sharpness, osc_0_0_0_0_0_p$spectral_sharpness),
            acf_regularity = c(osc_1_0_0_0_0_e$acf_regularity, osc_0_0_0_0_0_e$acf_regularity,
                               osc_1_0_0_0_0_p$acf_regularity, osc_0_0_0_0_0_p$acf_regularity)
          )
          
          # # ==== interpretation
          # # ==== acf_peak_val > 0.3–0.4 typically indicates clear oscillations (like your 39518 example)
          # # ==== acf_peak_lag gives you the approximate period
          # # ==== spec_conc > 0.1–0.15 suggests concentrated periodic power vs diffuse noise
          
          # Append to global results
          all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
        }
      }
      
      if(dim(all_comparison_results_reps)[1]>0){
        # ====== EPITHELIAL SCORE
        # JS Divergence comparisons
        # 1: Control vs tregs OFF (baseline effect of ROS/injury)
        all_comparison_results_reps$d_ctrl_vs_tregs_off_pathogen_e  = calculate_js_divergence(scores_1_0_0_0_0_e_keep, scores_0_0_0_0_0_e_keep, method = "fd")[1]
        
        # ===== PATHOGEN ABUNDANCE
        # JS Divergence comparisons
        # 1: Control vs tregs OFF (baseline effect of ROS/injury)
        all_comparison_results_reps$d_ctrl_vs_tregs_off_pathogen_p  = calculate_js_divergence(scores_1_0_0_0_0_p_keep, scores_0_0_0_0_0_p_keep, method = "fd")[1]
        
        # ====== Mean scores - epithelial score
        # Control
        all_comparison_results_reps$mean_ctrl_pathogen_e = mean(scores_1_0_0_0_0_e_keep)
        # Test - tregs OFF
        all_comparison_results_reps$mean_tregs_off_pathogen_e = mean(scores_0_0_0_0_0_e_keep)
        
        # ====== SD - epithelial score
        # Control
        all_comparison_results_reps$sd_ctrl_pathogen_e = sd(scores_1_0_0_0_0_e_keep)
        
        # Test - tregs OFF
        all_comparison_results_reps$sd_tregs_off_pathogen_e = sd(scores_0_0_0_0_0_e_keep)
        
        # ====== Mean scores - pathogen abundance
        # Control
        all_comparison_results_reps$mean_ctrl_pathogen_p = mean(scores_1_0_0_0_0_p_keep)
        
        # Test - tregs OFF
        all_comparison_results_reps$mean_tregs_off_pathogen_p = mean(scores_0_0_0_0_0_p_keep)
        
        # ====== SD - pathogen abundance
        # Control
        all_comparison_results_reps$sd_ctrl_pathogen_p = sd(scores_1_0_0_0_0_p_keep)
        
        # Test - tregs OFF
        all_comparison_results_reps$sd_tregs_off_pathogen_p = sd(scores_0_0_0_0_0_p_keep)
        
        all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
        processed_indices      = c(processed_indices, i)
      }
    }
    # Save after every 10 parameter sets (if total is > 10) or save all (if total is <= 10)
    if((length(inds2read) > 100 && i_idx %% 100 == 0) || i_idx==length(loop_over)){
      message("Saving intermediate results after ", i_idx, " parameter sets...")
      
      # Update the list of read indices
      if(!file.exists('./ids_read_abm_ros_vs_ctrl.rds')){
        saveRDS(processed_indices, './ids_read_abm_ros_vs_ctrl.rds')
      }else{
        inds_read_old = readRDS('./ids_read_abm_ros_vs_ctrl.rds')
        inds_read_updated = c(inds_read_old, processed_indices)
        saveRDS(inds_read_updated, './ids_read_abm_ros_vs_ctrl.rds')
      }
      
      # Update the results data
      if(!file.exists('./data_cpp_read_abm_ros_vs_ctrl.rds')){
        saveRDS(all_comparison_results, './data_cpp_read_abm_ros_vs_ctrl.rds')
      }else{
        all_comparison_results_old = readRDS('./data_cpp_read_abm_ros_vs_ctrl.rds')
        all_comparison_results_combined = rbind(all_comparison_results_old, all_comparison_results)
        saveRDS(all_comparison_results_combined, './data_cpp_read_abm_ros_vs_ctrl.rds')
      }
      
      # Reset for next batch
      all_comparison_results = data.frame()
      processed_indices = c()
      
      message("Intermediate save complete with 10 more param_ids.")
    }
  }
  
  # Final save for any remaining results
  message("Saving final results. Total param_sets processed: ", length(inds2read))
  
  if(!file.exists('./ids_read_abm_ros_vs_ctrl.rds')){
    saveRDS(processed_indices, './ids_read_abm_ros_vs_ctrl.rds')
  }else{
    inds_read_old = readRDS('./ids_read_abm_ros_vs_ctrl.rds')
    inds_read_final = c(inds_read_old, processed_indices)
    saveRDS(inds_read_final, './ids_read_abm_ros_vs_ctrl.rds')
  }
  
  if(!file.exists('./data_cpp_read_abm_ros_vs_ctrl.rds')){
    saveRDS(all_comparison_results, './data_cpp_read_abm_ros_vs_ctrl.rds')
  }else{
    all_comparison_results_old = readRDS('./data_cpp_read_abm_ros_vs_ctrl.rds')
    all_comparison_results_final = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results_final, './data_cpp_read_abm_ros_vs_ctrl.rds')
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




