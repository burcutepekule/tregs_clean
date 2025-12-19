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

# Control scenario (reference)
# control_1_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0
files_ctrl = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_control_1_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds$", full.names = TRUE)

# Test scenarios with different ros_level and pat_level combinations
# We'll look for all combinations of ros_level (0,1) and pat_level (1,2,3) with control=0, tregs=0
test_scenarios = expand.grid(
  ros_level = c(0, 1),
  pat_level = c(1, 2, 3)
)

# Get all test files
test_files_list = list()
for (i in 1:nrow(test_scenarios)) {
  ros_lvl = test_scenarios$ros_level[i]
  pat_lvl = test_scenarios$pat_level[i]
  pattern = sprintf("^longitudinal_df_param_set_id_\\d+\\_control_0_sterile_0_macspec_0_tregs_0_ros_level_%d_pat_level_%d_trnd_0.rds$",
                    ros_lvl, pat_lvl)
  files = list.files(path, pattern = pattern, full.names = TRUE)
  test_files_list[[paste0("ros_", ros_lvl, "_pat_", pat_lvl)]] = files
}

# Extract indices from control files
indices_ctrl = str_extract(basename(files_ctrl), "param_set_id_(\\d+)_", group = 1) |> as.numeric()

# Extract indices from all test scenarios and find common indices
all_test_indices = list()
for (scenario_name in names(test_files_list)) {
  files = test_files_list[[scenario_name]]
  if (length(files) > 0) {
    all_test_indices[[scenario_name]] = str_extract(basename(files), "param_set_id_(\\d+)_", group = 1) |> as.numeric()
  }
}

# Find common indices across control and all test scenarios
indices = Reduce(intersect, c(list(indices_ctrl), all_test_indices))

message("Found ", length(indices), " common parameter sets across all scenarios")

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('./ids_read_abm_ros_vs_ctrl_patros.rds')){
  inds_read = readRDS('./ids_read_abm_ros_vs_ctrl_patros.rds')
}else{
  inds_read = c()
}

inds2read = sort(setdiff(indices,inds_read))
message("New parameter sets to read: ", length(inds2read))

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

    # Check file sizes for control and all test scenarios
    files_to_check = c(
      paste0(path, 'longitudinal_df_param_set_id_', i, '_control_1_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds')
    )

    # Add all test scenario files to check
    for (row_idx in 1:nrow(test_scenarios)) {
      ros_lvl = test_scenarios$ros_level[row_idx]
      pat_lvl = test_scenarios$pat_level[row_idx]
      files_to_check = c(files_to_check,
                        paste0(path, 'longitudinal_df_param_set_id_', i,
                              '_control_0_sterile_0_macspec_0_tregs_0_ros_level_', ros_lvl,
                              '_pat_level_', pat_lvl, '_trnd_0.rds'))
    }

    if(any(!file.exists(files_to_check))){
      message("Some files missing for param_set_", i, " - skipping")
      processed_indices = c(processed_indices, i) # add and skip
    }else if(any(file.info(files_to_check)$size<100000)){
      processed_indices = c(processed_indices, i) # add and skip
      message("Skipped one due to small file size")
    }else{
      # Read control file
      results_ctrl = readRDS(paste0(path, 'longitudinal_df_param_set_id_', i,
                                   '_control_1_sterile_0_macspec_0_tregs_0_ros_level_0_pat_level_1_trnd_0.rds'))

      # Read all test scenario files
      test_results_list = list()
      for (row_idx in 1:nrow(test_scenarios)) {
        ros_lvl = test_scenarios$ros_level[row_idx]
        pat_lvl = test_scenarios$pat_level[row_idx]
        scenario_name = paste0("ros_", ros_lvl, "_pat_", pat_lvl)
        test_results_list[[scenario_name]] = readRDS(
          paste0(path, 'longitudinal_df_param_set_id_', i,
                '_control_0_sterile_0_macspec_0_tregs_0_ros_level_', ros_lvl,
                '_pat_level_', pat_lvl, '_trnd_0.rds')
        )
      }

      # Combine all results for analysis
      results = do.call(rbind, c(list(results_ctrl), test_results_list))

      full_data_comparison = results %>%
        dplyr::select(param_set_id, control, sterile, macspec_on, tregs_on, randomize_tregs,
                     ros_level, pat_level, rep_id, t, time_ss, epithelial_score, pathogen)

      min_reps  = min(full_data_comparison$rep_id)
      max_reps  = max(full_data_comparison$rep_id)
      t_max_ind = max(full_data_comparison$t)

      all_comparison_results_reps = data.frame()

      for (rep in min_reps:max_reps) {

        # Control scenario
        full_data_ctrl = full_data_comparison %>%
          dplyr::filter(rep_id==rep & control==1 & sterile==0 & macspec_on==0 &
                       tregs_on==0 & randomize_tregs==0 & ros_level==0 & pat_level==1)

        # Process each test scenario
        for (row_idx in 1:nrow(test_scenarios)) {
          ros_lvl = test_scenarios$ros_level[row_idx]
          pat_lvl = test_scenarios$pat_level[row_idx]

          full_data_test = full_data_comparison %>%
            dplyr::filter(rep_id==rep & control==0 & sterile==0 & macspec_on==0 &
                         tregs_on==0 & randomize_tregs==0 &
                         ros_level==ros_lvl & pat_level==pat_lvl)

          # --- Steady-state detection ---
          time_ss_ctrl_e = as.numeric(steady_state_idx(full_data_ctrl$epithelial_score))
          time_ss_ctrl_p = as.numeric(steady_state_idx(full_data_ctrl$pathogen))

          time_ss_test_e = as.numeric(steady_state_idx(full_data_test$epithelial_score))
          time_ss_test_p = as.numeric(steady_state_idx(full_data_test$pathogen))

          time_ss_vec = c(time_ss_ctrl_e, time_ss_ctrl_p, time_ss_test_e, time_ss_test_p)

          if(!any(is.na(time_ss_vec))){

            # Extract steady-state scores
            scores_ctrl_p = full_data_ctrl$pathogen[time_ss_ctrl_p:t_max_ind]
            scores_test_p = full_data_test$pathogen[time_ss_test_p:t_max_ind]

            scores_ctrl_e = full_data_ctrl$epithelial_score[time_ss_ctrl_e:t_max_ind]
            scores_test_e = full_data_test$epithelial_score[time_ss_test_e:t_max_ind]

            # Compute oscillation metrics for each signal
            osc_ctrl_e = compute_oscillation_metrics(scores_ctrl_e)
            osc_test_e = compute_oscillation_metrics(scores_test_e)
            osc_ctrl_p = compute_oscillation_metrics(scores_ctrl_p)
            osc_test_p = compute_oscillation_metrics(scores_test_p)

            # --- Tabulate comparison results ---
            comparison_results = data.frame(
              param_set_id = i,
              replicate_id = rep,
              ros_level    = c(0, ros_lvl, 0, ros_lvl),
              pat_level    = c(1, pat_lvl, 1, pat_lvl),
              control      = c(1, 0, 1, 0),
              injury_type  = c("pathogenic", "pathogenic", "pathogenic", "pathogenic"),
              score_type   = c("epithelium", "epithelium", "pathogen", "pathogen"),
              macspec_on   = c(0, 0, 0, 0),
              tregs_on     = c(0, 0, 0, 0),
              tregs_rnd    = c(0, 0, 0, 0),
              ss_start     = c(time_ss_ctrl_e, time_ss_test_e, time_ss_ctrl_p, time_ss_test_p),
              mean_score   = c(mean(scores_ctrl_e), mean(scores_test_e),
                              mean(scores_ctrl_p), mean(scores_test_p)),
              sd_score     = c(sd(scores_ctrl_e), sd(scores_test_e),
                              sd(scores_ctrl_p), sd(scores_test_p)),
              # Oscillation metrics
              acf_peak_lag = c(osc_ctrl_e$acf_peak_lag, osc_test_e$acf_peak_lag,
                              osc_ctrl_p$acf_peak_lag, osc_test_p$acf_peak_lag),
              acf_peak_val = c(osc_ctrl_e$acf_peak_value, osc_test_e$acf_peak_value,
                              osc_ctrl_p$acf_peak_value, osc_test_p$acf_peak_value),
              spec_conc    = c(osc_ctrl_e$spectral_concentration, osc_test_e$spectral_concentration,
                              osc_ctrl_p$spectral_concentration, osc_test_p$spectral_concentration),
              spec_sharp   = c(osc_ctrl_e$spectral_sharpness, osc_test_e$spectral_sharpness,
                              osc_ctrl_p$spectral_sharpness, osc_test_p$spectral_sharpness),
              acf_regularity = c(osc_ctrl_e$acf_regularity, osc_test_e$acf_regularity,
                                osc_ctrl_p$acf_regularity, osc_test_p$acf_regularity)
            )

            # Compute JS divergences for this ros_level and pat_level combination
            comparison_results$d_ctrl_vs_test_e = calculate_js_divergence(scores_ctrl_e, scores_test_e, method = "fd")[1]
            comparison_results$d_ctrl_vs_test_p = calculate_js_divergence(scores_ctrl_p, scores_test_p, method = "fd")[1]

            # Append to results
            all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
          }
        }
      }

      if(dim(all_comparison_results_reps)[1]>0){
        all_comparison_results = bind_rows(all_comparison_results, all_comparison_results_reps)
        processed_indices      = c(processed_indices, i)
      }
    }

    # Save after every 10 parameter sets (if total is > 10) or save all (if total is <= 10)
    if((length(inds2read) > (10*n2) && i_idx %% (10*n2) == 0) || i_idx==length(loop_over)){ #10*n2 to break sync between multiple cores

      message("Saving intermediate results after ", i_idx, " parameter sets...")

      # Update the list of read indices
      if(!file.exists('./ids_read_abm_ros_vs_ctrl_patros.rds')){
        saveRDS(processed_indices, './ids_read_abm_ros_vs_ctrl_patros.rds')
      }else{
        inds_read_old = readRDS('./ids_read_abm_ros_vs_ctrl_patros.rds')
        inds_read_updated = c(inds_read_old, processed_indices)
        saveRDS(inds_read_updated, './ids_read_abm_ros_vs_ctrl_patros.rds')
      }

      # Update the results data
      if(!file.exists('./data_cpp_read_abm_ros_vs_ctrl_patros.rds')){
        saveRDS(all_comparison_results, './data_cpp_read_abm_ros_vs_ctrl_patros.rds')
      }else{
        all_comparison_results_old = readRDS('./data_cpp_read_abm_ros_vs_ctrl_patros.rds')
        all_comparison_results_combined = rbind(all_comparison_results_old, all_comparison_results)
        saveRDS(all_comparison_results_combined, './data_cpp_read_abm_ros_vs_ctrl_patros.rds')
      }

      # Reset for next batch
      all_comparison_results = data.frame()
      processed_indices = c()

      message(paste0("Intermediate save complete with ",10*n2," more param_ids."))
    }
  }

  # Final save for any remaining results
  message("Saving final results. Total param_sets processed: ", length(inds2read))

  if(!file.exists('./ids_read_abm_ros_vs_ctrl_patros.rds')){
    saveRDS(processed_indices, './ids_read_abm_ros_vs_ctrl_patros.rds')
  }else{
    inds_read_old = readRDS('./ids_read_abm_ros_vs_ctrl_patros.rds')
    inds_read_final = c(inds_read_old, processed_indices)
    saveRDS(inds_read_final, './ids_read_abm_ros_vs_ctrl_patros.rds')
  }

  if(!file.exists('./data_cpp_read_abm_ros_vs_ctrl_patros.rds')){
    saveRDS(all_comparison_results, './data_cpp_read_abm_ros_vs_ctrl_patros.rds')
  }else{
    all_comparison_results_old = readRDS('./data_cpp_read_abm_ros_vs_ctrl_patros.rds')
    all_comparison_results_final = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results_final, './data_cpp_read_abm_ros_vs_ctrl_patros.rds')
  }

}else{
  message("No new pts added.")
}



