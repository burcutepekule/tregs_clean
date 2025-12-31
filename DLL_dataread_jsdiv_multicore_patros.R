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

# Define the ros and pat value ranges
ros_vals = c(0, 0.1, 0.25, seq(0.5,10,0.5))
pat_vals = c(1:15)

path = "/Users/burcutepekule/Desktop/DLL/sim_abm/"

# File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
# Dynamically create file lists for all ros/pat combinations
all_files = list()
all_indices = list()

for (ros in ros_vals) {
  for (pat in pat_vals) {
    var_name = paste0("files_", ros, "_", pat)
    pattern_str = paste0("^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_0_ros_level_", ros, "_pat_level_", pat, "_trnd_0.rds$")
    files_temp = list.files(path, pattern = pattern_str, full.names = TRUE)
    all_files[[var_name]] = files_temp

    # Extract indices
    indices_name = paste0("indices_", ros, "_", pat)
    all_indices[[indices_name]] = str_extract(basename(files_temp), "\\d+") |> as.numeric()
  }
}

# Find common indices across all combinations
indices = Reduce(intersect, all_indices)

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

    # Dynamically create file paths to check
    files_to_check = c()
    for (ros in ros_vals) {
      for (pat in pat_vals) {
        files_to_check = c(files_to_check,
          paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_', ros, '_pat_level_', pat, '_trnd_0.rds'))
      }
    }

    if(any(file.info(files_to_check)$size<10000)){
      processed_indices = c(processed_indices, i) #add and skip
      message("Skipped one due to file size.")
    }else{

      # Dynamically read all RDS files for this parameter set
      results_list = list()
      for (ros in ros_vals) {
        for (pat in pat_vals) {
          var_name = paste0("results_", ros, "_", pat)
          file_path = paste0(path, 'longitudinal_df_param_set_id_', i, '_sterile_0_macspec_0_tregs_0_ros_level_', ros, '_pat_level_', pat, '_trnd_0.rds')
          results_list[[var_name]] = readRDS(file_path)
        }
      }

      # Combine all results
      results = do.call(rbind, results_list)
      
      full_data_comparison = results %>% dplyr::select(param_set_id, sterile, macspec_on, tregs_on, 
                                                       randomize_tregs, ros_level, pat_level, rep_id, 
                                                       t, time_ss_e, time_ss_p, epithelial_score, pathogen)
      min_reps  = min(full_data_comparison$rep_id)
      max_reps  = max(full_data_comparison$rep_id)
      t_max_ind = max(full_data_comparison$t)

      # Dynamically initialize score keeping variables for all combinations
      # PATHOGEN
      for (ros in ros_vals) {
        for (pat in pat_vals) {
          assign(paste0("scores_", ros, "_", pat, "_p_keep"), c())
        }
      }

      # EPITHELIUM
      for (ros in ros_vals) {
        for (pat in pat_vals) {
          assign(paste0("scores_", ros, "_", pat, "_e_keep"), c())
        }
      }
      
      all_comparison_results_reps = data.frame()
      
      for (rep in min_reps:max_reps) {

        # Dynamically filter data and compute steady states for all combinations
        time_ss_vec = c()
        for (ros in ros_vals) {
          for (pat in pat_vals) {
            # Filter data
            var_name = paste0("full_data_comparison_scores_", ros, "_", pat)
            assign(var_name, full_data_comparison %>%
              dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==0 &
                           ros_level==ros & pat_level==pat & randomize_tregs==0))

            # Compute steady state for epithelial score
            time_ss_e_var = paste0("time_ss_", ros, "_", pat, "_e")
            # time_ss_e_val = as.numeric(steady_state_idx(get(var_name)$epithelial_score))
            time_ss_e_val = as.numeric(unique(get(var_name)$time_ss_e))
            
            assign(time_ss_e_var, time_ss_e_val)
            time_ss_vec = c(time_ss_vec, time_ss_e_val)

            # Compute steady state for pathogen
            time_ss_p_var = paste0("time_ss_", ros, "_", pat, "_p")
            # time_ss_p_val = as.numeric(steady_state_idx(get(var_name)$pathogen))
            time_ss_p_val = as.numeric(unique(get(var_name)$time_ss_p))
            
            assign(time_ss_p_var, time_ss_p_val)
            time_ss_vec = c(time_ss_vec, time_ss_p_val)
          }
        }
        
        if(!any(is.na(time_ss_vec))){

          # Dynamically extract scores and accumulate them
          for (ros in ros_vals) {
            for (pat in pat_vals) {
              # Get the filtered data and steady state times
              scores_df = get(paste0("full_data_comparison_scores_", ros, "_", pat))
              time_ss_p = get(paste0("time_ss_", ros, "_", pat, "_p"))
              time_ss_e = get(paste0("time_ss_", ros, "_", pat, "_e"))

              # Extract pathogen scores
              scores_p_var = paste0("scores_", ros, "_", pat, "_p")
              scores_p_val = scores_df$pathogen[time_ss_p:t_max_ind]
              assign(scores_p_var, scores_p_val)

              # Accumulate pathogen scores
              scores_p_keep_var = paste0("scores_", ros, "_", pat, "_p_keep")
              assign(scores_p_keep_var, c(get(scores_p_keep_var), scores_p_val))

              # Extract epithelial scores
              scores_e_var = paste0("scores_", ros, "_", pat, "_e")
              scores_e_val = scores_df$epithelial_score[time_ss_e:t_max_ind]
              assign(scores_e_var, scores_e_val)

              # Accumulate epithelial scores
              scores_e_keep_var = paste0("scores_", ros, "_", pat, "_e_keep")
              assign(scores_e_keep_var, c(get(scores_e_keep_var), scores_e_val))
            }
          }
          
          # Compute oscillation metrics for each signal 
          # for (ros in ros_vals) {
          #   for (pat in pat_vals) {
          #     osc_e_var = paste0("osc_", ros, "_", pat, "_e")
          #     assign(osc_e_var, compute_oscillation_metrics(get(paste0("scores_", ros, "_", pat, "_e"))))
          # 
          #     osc_p_var = paste0("osc_", ros, "_", pat, "_p")
          #     assign(osc_p_var, compute_oscillation_metrics(get(paste0("scores_", ros, "_", pat, "_p"))))
          #   }
          # }
          
          # --- Tabulate all comparisons ---
          # Dynamically build the comparison_results data frame
          n_combinations = length(ros_vals) * length(pat_vals)
          n_total_rows = 2 * n_combinations  # epithelium + pathogen

          # Build vectors for all combinations
          ros_vec_e = rep(ros_vals, each = length(pat_vals))
          pat_vec_e = rep(pat_vals, times = length(ros_vals))
          ros_vec_p = ros_vec_e
          pat_vec_p = pat_vec_e

          # Combine for epithelium and pathogen
          ros_vec_all = c(ros_vec_e, ros_vec_p)
          pat_vec_all = c(pat_vec_e, pat_vec_p)

          # Build ss_start, mean_score, and sd_score vectors dynamically
          ss_start_vec = c()
          mean_score_vec = c()
          sd_score_vec = c()

          # Epithelium first
          for (ros in ros_vals) {
            for (pat in pat_vals) {
              ss_start_vec = c(ss_start_vec, get(paste0("time_ss_", ros, "_", pat, "_e")))
              mean_score_vec = c(mean_score_vec, mean(get(paste0("scores_", ros, "_", pat, "_e"))))
              sd_score_vec = c(sd_score_vec, sd(get(paste0("scores_", ros, "_", pat, "_e"))))
            }
          }

          # Pathogen second
          for (ros in ros_vals) {
            for (pat in pat_vals) {
              ss_start_vec = c(ss_start_vec, get(paste0("time_ss_", ros, "_", pat, "_p")))
              mean_score_vec = c(mean_score_vec, mean(get(paste0("scores_", ros, "_", pat, "_p"))))
              sd_score_vec = c(sd_score_vec, sd(get(paste0("scores_", ros, "_", pat, "_p"))))
            }
          }

          comparison_results = data.frame(
            param_set_id = i,
            replicate_id = rep,
            injury_type  = rep("pathogenic", n_total_rows),
            score_type   = c(rep("epithelium", n_combinations), rep("pathogen", n_combinations)),
            macspec_on   = rep(0, n_total_rows),
            tregs_on     = rep(0, n_total_rows),
            tregs_rnd    = rep(0, n_total_rows),
            ros_level    = ros_vec_all,
            pat_level    = pat_vec_all,
            ss_start     = ss_start_vec,
            mean_score   = mean_score_vec,
            sd_score     = sd_score_vec
          )
          
          # # ==== interpretation
          # # ==== acf_peak_val > 0.3–0.4 typically indicates clear oscillations (like your 39518 example)
          # # ==== acf_peak_lag gives you the approximate period
          # # ==== spec_conc > 0.1–0.15 suggests concentrated periodic power vs diffuse noise
          
          # Append to global results
          all_comparison_results_reps = bind_rows(all_comparison_results_reps, comparison_results)
        }
        else{
          # print(time_ss_vec)
          message("Skipping replicate ", rep, " for param_set_", i, " due to NA in time_ss_vec")
          processed_indices = c(processed_indices, i)
          next
        }
      }
      
      if(dim(all_comparison_results_reps)[1]>0){

        # JS Divergence - ALL POSSIBLE COMBINATIONS
        # All comparisons between all conditions with explicit naming

        # # ==== EPITHELIAL SCORE
        # for (r1 in ros_vals) {
        #   for (p1 in pat_vals) {
        #     for (r2 in ros_vals) {
        #       for (p2 in pat_vals) {
        #         if (r1 != r2 || p1 != p2) {
        #           col_name = paste0("d_ros", r1, "_pat", p1, "_treg0_vs_ros", r2, "_pat", p2, "_treg0_e")
        #           scores1 = get(paste0("scores_", r1, "_", p1, "_e_keep"))
        #           scores2 = get(paste0("scores_", r2, "_", p2, "_e_keep"))
        #           all_comparison_results_reps[[col_name]] = calculate_js_divergence(scores1, scores2, method = "fd")[1]
        #         }
        #       }
        #     }
        #   }
        # }
        # 
        # # ==== PATHOGEN ABUNDANCE
        # for (r1 in ros_vals) {
        #   for (p1 in pat_vals) {
        #     for (r2 in ros_vals) {
        #       for (p2 in pat_vals) {
        #         if (r1 != r2 || p1 != p2) {
        #           col_name = paste0("d_ros", r1, "_pat", p1, "_treg0_vs_ros", r2, "_pat", p2, "_treg0_p")
        #           scores1 = get(paste0("scores_", r1, "_", p1, "_p_keep"))
        #           scores2 = get(paste0("scores_", r2, "_", p2, "_p_keep"))
        #           all_comparison_results_reps[[col_name]] = calculate_js_divergence(scores1, scores2, method = "fd")[1]
        #         }
        #       }
        #     }
        #   }
        # }
        
        # ====== Mean scores - epithelial score (dynamically)
        for (ros in ros_vals) {
          for (pat in pat_vals) {
            col_name = paste0("mean_ros", ros, "_pat", pat, "_treg0_e")
            all_comparison_results_reps[[col_name]] = mean(get(paste0("scores_", ros, "_", pat, "_e_keep")))
          }
        }

        # ====== Mean scores - pathogen abundance (dynamically)
        for (ros in ros_vals) {
          for (pat in pat_vals) {
            col_name = paste0("mean_ros", ros, "_pat", pat, "_treg0_p")
            all_comparison_results_reps[[col_name]] = mean(get(paste0("scores_", ros, "_", pat, "_p_keep")))
          }
        }
        
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



