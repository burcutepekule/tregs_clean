is_under_control = function(epithelial_health, pathogen_load) {
  return(epithelial_health > 150*0.75 & pathogen_load < 10)
}

message("Re-processing param_set_", param_id)

# Dynamically create file paths to check
files_to_check = c()
for (ros in ros_vals) {
  for (pat in pat_vals) {
    files_to_check = c(files_to_check,
                       paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_',sterile_in,'_macspec_0_tregs_',tregs_on_in,'_ros_level_', ros, '_pat_level_', pat, '_trnd_0_optidx_', opt_idx,'.rds'))
  }
}

files2read = which(file.info(files_to_check)$size>10000)

# Dynamically read all RDS files for this parameter set
results_list = list()
for (ros in ros_vals) {
  for (pat in pat_vals) {
    # print(paste0('Reading ros ',ros,' pat ',pat))
    var_name = paste0("results_", ros, "_", pat)
    file_path = paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_',sterile_in,'_macspec_0_tregs_',tregs_on_in,'_ros_level_', ros, '_pat_level_', pat, '_trnd_0_optidx_', opt_idx,'.rds')
    if(file.exists(file_path) & file.info(file_path)$size>10000){
      results_list[[var_name]] = readRDS(file_path)
    }
  }
}

# Combine all results
results = do.call(rbind, results_list)

full_data_comparison = results %>% dplyr::select(param_set_id, sterile, macspec_on, tregs_on, 
                                                 randomize_tregs, ros_level, pat_level, rep_id, 
                                                 t, time_ss_e, time_ss_p, epithelial_score, pathogen)
min_reps  = min(full_data_comparison$rep_id)
max_reps  = min(max(rep_ind_vec),max(full_data_comparison$rep_id))
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
all_comparison_results      = data.frame()

for (rep in min_reps:max_reps) {
  print(paste0('Processing rep ',rep))
  # Dynamically filter data and compute steady states for all combinations
  time_ss_vec = c()
  for (ros in ros_vals) {
    for (pat in pat_vals) {
      # Filter data
      var_name = paste0("full_data_comparison_scores_", ros, "_", pat)
      assign(var_name, full_data_comparison %>%
               dplyr::filter(rep_id==rep & sterile==0 & macspec_on==0 & tregs_on==tregs_on_in &
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
  
  # if(!any(is.na(time_ss_vec))){
  if(length(time_ss_vec)>0){
    
    # Dynamically extract scores and accumulate them
    for (ros in ros_vals) {
      for (pat in pat_vals) {
        # Get the filtered data and steady state times
        scores_df = get(paste0("full_data_comparison_scores_", ros, "_", pat))
        if(dim(scores_df)[1]>0){
          time_ss_p = get(paste0("time_ss_", ros, "_", pat, "_p"))
          time_ss_e = get(paste0("time_ss_", ros, "_", pat, "_e"))
          
          if(is.na(time_ss_p)){time_ss_p = t_max_ind-50}
          if(is.na(time_ss_e)){time_ss_e = t_max_ind-50}
          
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
    

    # # Build vectors for all combinations
    # ros_vec_e = rep(ros_vals, each = length(pat_vals))
    # pat_vec_e = rep(pat_vals, times = length(ros_vals))
    # ros_vec_p = ros_vec_e
    # pat_vec_p = pat_vec_e
    # 
    # # Combine for epithelium and pathogen
    # ros_vec_all = c(ros_vec_e, ros_vec_p)
    # pat_vec_all = c(pat_vec_e, pat_vec_p)
    
    # Build ss_start, mean_score, and sd_score vectors dynamically
    ss_start_vec = c()
    mean_score_vec = c()
    sd_score_vec = c()
    
    ros_vec_all = c()
    pat_vec_all = c()
    # Epithelium first
    for (ros in ros_vals) {
      for (pat in pat_vals) {
        var_time <- paste0("time_ss_", ros, "_", pat, "_e")
        var_scores <- paste0("scores_", ros, "_", pat, "_e")
        
        if (exists(var_time) && exists(var_scores)) {
          ss_start_vec <- c(ss_start_vec, get(var_time))
          mean_score_vec <- c(mean_score_vec, mean(get(var_scores)))
          sd_score_vec <- c(sd_score_vec, sd(get(var_scores)))
          ros_vec_all = c(ros_vec_all, ros)
          pat_vec_all = c(pat_vec_all, pat)
        }
      }
    }
    
    # Pathogen second
    for (ros in ros_vals) {
      for (pat in pat_vals) {
        var_time <- paste0("time_ss_", ros, "_", pat, "_p")
        var_scores <- paste0("scores_", ros, "_", pat, "_p")
        
        if (exists(var_time) && exists(var_scores)) {
          ss_start_vec = c(ss_start_vec, get(var_time))
          mean_score_vec = c(mean_score_vec, mean(get(var_scores)))
          sd_score_vec = c(sd_score_vec, sd(get(var_scores)))
        }
      }
    }
    
    # --- Tabulate all comparisons ---
    # Dynamically build the comparison_results data frame
    n_total_rows   = length(ss_start_vec)  # epithelium + pathogen
    n_combinations = n_total_rows/2
    
    comparison_results = data.frame(
      param_set_id = param_id,
      replicate_id = rep,
      injury_type  = rep("pathogenic", n_total_rows),
      score_type   = c(rep("epithelium", n_combinations), rep("pathogen", n_combinations)),
      macspec_on   = rep(0, n_total_rows),
      tregs_on     = rep(tregs_on_in, n_total_rows),
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
    message("Skipping replicate ", rep, " for param_set_", param_id, " due to NA in time_ss_vec")
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
}



df_steps = all_comparison_results %>% dplyr::filter(injury_type == 'pathogenic')

# Dynamically create the controlled status columns for all ros/pat combinations
for (ros in ros_vals) {
  for (pat in pat_vals) {
    col_name = paste0("ros", ros, "_pat", pat, "_controlled")
    mean_e_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_e")
    mean_p_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_p")
    
    df_steps = df_steps %>%
      dplyr::mutate(!!col_name := is_under_control(.data[[mean_e_col]], .data[[mean_p_col]]))
  }
}

row = df_steps[df_steps$param_set_id == param_id, ]
row_vec = unlist(row)

# Build the control matrix
control_matrix = matrix(NA, nrow = length(pat_vals), ncol = length(ros_vals))
rownames(control_matrix) = paste0("pat", pat_vals)
colnames(control_matrix) = paste0("ros", ros_vals)

for (i in seq_along(pat_vals)) {
  for (j in seq_along(ros_vals)) {
    pat = pat_vals[i]
    ros = ros_vals[j]
    col_name = paste0("ros", ros, "_pat", pat, "_controlled")
    control_matrix[i, j] = as.logical(row[[col_name]][1])
  }
}


