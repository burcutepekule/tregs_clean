# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
# Update injury site
injury_site_updated = which(epithelium$level_injury > 0)

# ========================================================================
# UPDATE SAMPs (from activated Tregs)
# ========================================================================
active_tregs = which(treg_phenotype == 1)
if (length(active_tregs) > 0) {
  # C++ version (20-50x faster)
  SAMPs = update_SAMPs_batch_cpp(
    SAMPs, active_tregs, treg_x, treg_y,
    treg_activity_SAMPs_binary, add_SAMPs, allow_tregs
  )
}

# ========================================================================
# UPDATE ROS (from M1 phagocytes)
# ========================================================================
M1_phagocytes = which(phagocyte_phenotype == 1)
if (length(M1_phagocytes) > 0) {
  coords = cbind(phagocyte_y[M1_phagocytes], phagocyte_x[M1_phagocytes])
  ROS[coords] = ROS[coords] + phagocyte_activity_ROS[M1_phagocytes]*add_ROS
}

# ========================================================================
# MOVE MICROBES
# ========================================================================
if (nrow(pathogen_coords) > 0) {
  dy = ifelse(pathogen_coords[, "y"] == 1,
              sample(c(1), size = nrow(pathogen_coords), replace = TRUE),
              sample(c(-1, 0, 1), size = nrow(pathogen_coords), replace = TRUE))
  dx = iszero_coordinates(dy)
  pathogen_coords[, "x"] = pmin(pmax(pathogen_coords[, "x"] + dx, 1), grid_size)
  pathogen_coords[, "y"] = pmin(pmax(pathogen_coords[, "y"] + dy, 1), grid_size)
}

if (nrow(commensal_coords) > 0) {
  dy = ifelse(commensal_coords[, "y"] == 1,
              sample(c(1), size = nrow(commensal_coords), replace = TRUE),
              sample(c(-1, 0, 1), size = nrow(commensal_coords), replace = TRUE))
  dx = iszero_coordinates(dy)
  commensal_coords[, "x"] = pmin(pmax(commensal_coords[, "x"] + dx, 1), grid_size)
  commensal_coords[, "y"] = pmin(pmax(commensal_coords[, "y"] + dy, 1), grid_size)
}

# ========================================================================
# PRE-CALCULATE MICROBE COUNTS AT EPITHELIUM (y=1)
# ========================================================================
pathogen_epithelium_counts = rep(0, grid_size)
if (nrow(pathogen_coords) > 0) {
  epithelium_pathogens = pathogen_coords[pathogen_coords[, "y"] == 1, , drop = FALSE]
  if (nrow(epithelium_pathogens) > 0) {
    pathogen_epithelium_counts = tabulate(epithelium_pathogens[, "x"], nbins = grid_size)
  }
}

commensal_epithelium_counts = rep(0, grid_size)
if (nrow(commensal_coords) > 0) {
  epithelium_commensals = commensal_coords[commensal_coords[, "y"] == 1, , drop = FALSE]
  if (nrow(epithelium_commensals) > 0) {
    commensal_epithelium_counts = tabulate(epithelium_commensals[, "x"], nbins = grid_size)
  }
}

# ========================================================================
# UPDATE DAMPs 
# ========================================================================
# DAMPs from epithelial injury
DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury*add_DAMPs

# DAMPs from microbes touching epithelium
DAMPs[1, ] = DAMPs[1, ] + logistic_scaled_0_to_5_quantized(
  pathogen_epithelium_counts + commensal_epithelium_counts
)*add_DAMPs

# ========================================================================
# UPDATE PAMPs 
# ========================================================================
# PAMPs from pathogen locations - toxins released by pathogens
PAMPs_add = matrix(0, nrow = nrow(PAMPs), ncol = ncol(PAMPs))
if (nrow(pathogen_coords) > 0) {
  pat_counts_tab = table(pathogen_coords[, "x"], pathogen_coords[, "y"])
  if (length(pat_counts_tab) > 0) {
    pat_x = as.numeric(rownames(pat_counts_tab))
    pat_y = as.numeric(colnames(pat_counts_tab))
    for (xi in seq_along(pat_x)) {
      for (yi in seq_along(pat_y)) {
        if (pat_counts_tab[xi, yi] > 0) {
          PAMPs_add[pat_y[yi], pat_x[xi]] = add_PAMPs *
            logistic_scaled_0_to_5_quantized(pat_counts_tab[xi, yi])
        }
      }
    }
  }
}
PAMPs = PAMPs + PAMPs_add

# ========================================================================
# DIFFUSE & DECAY SIGNALS
# C++ ACCELERATION: Matrix diffusion operations
# ========================================================================
# C++ version (5-10x faster)
DAMPs = diffuse_matrix_cpp(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs, reflect_top=FALSE)
SAMPs = diffuse_matrix_cpp(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs, reflect_top=FALSE)
PAMPs = diffuse_matrix_cpp(PAMPs, diffusion_speed_PAMPs, max_cell_value_PAMPs, reflect_top=FALSE) 
ROS   = diffuse_matrix_cpp(ROS, diffusion_speed_ROS, max_cell_value_ROS, reflect_top=FALSE)

DAMPs = DAMPs - DAMPs_decay*DAMPs
SAMPs = SAMPs - SAMPs_decay*SAMPs
PAMPs = PAMPs - PAMPs_decay*PAMPs
ROS   = ROS - ros_decay*ROS

# ========================================================================
# RECRUIT MACROPHAGES FROM BORDERS (based on danger signal)
# ========================================================================
if (recruitment_rate_danger > 0 && length(phagocyte_x) < max_total_phagocytes) {
  danger_signal_grid = DAMPs + PAMPs
  
  # Recruit from BOTTOM border (y = grid_size)
  for (x in 1:grid_size) {
    danger_at_border = danger_signal_grid[grid_size, x]
    n_recruit = rpois(1, lambda = recruitment_rate_danger*danger_at_border)
    
    if (n_recruit > 0) {
      phagocyte_x = c(phagocyte_x, rep(x, n_recruit))
      phagocyte_y = c(phagocyte_y, rep(grid_size, n_recruit))
      phagocyte_pathogens_engulfed = c(phagocyte_pathogens_engulfed, rep(0, n_recruit))
      phagocyte_commensals_engulfed = c(phagocyte_commensals_engulfed, rep(0, n_recruit))
      phagocyte_num_times_activated = c(phagocyte_num_times_activated, rep(0, n_recruit))
      phagocyte_phenotype = c(phagocyte_phenotype, rep(0, n_recruit))  # Start as M0
      phagocyte_activity_ROS = c(phagocyte_activity_ROS, rep(activity_ROS_M0_baseline, n_recruit))
      phagocyte_activity_engulf = c(phagocyte_activity_engulf, rep(activity_engulf_M0_baseline, n_recruit))
      phagocyte_active_age = c(phagocyte_active_age, rep(0, n_recruit))
      
      # Extend bacteria registry matrix
      new_registry = matrix(0, nrow = n_recruit, ncol = cc_phagocyte)
      phagocyte_bacteria_registry = rbind(phagocyte_bacteria_registry, new_registry)
    }
  }
  
  # Recruit from LEFT border (x = 1, y > 1 to avoid epithelium)
  for (y in 2:grid_size) {
    danger_at_border = danger_signal_grid[y, 1]
    n_recruit = rpois(1, lambda = recruitment_rate_danger*danger_at_border)
    
    if (n_recruit > 0) {
      phagocyte_x = c(phagocyte_x, rep(1, n_recruit))
      phagocyte_y = c(phagocyte_y, rep(y, n_recruit))
      phagocyte_pathogens_engulfed = c(phagocyte_pathogens_engulfed, rep(0, n_recruit))
      phagocyte_commensals_engulfed = c(phagocyte_commensals_engulfed, rep(0, n_recruit))
      phagocyte_num_times_activated = c(phagocyte_num_times_activated, rep(0, n_recruit))
      phagocyte_phenotype = c(phagocyte_phenotype, rep(0, n_recruit))
      phagocyte_activity_ROS = c(phagocyte_activity_ROS, rep(activity_ROS_M0_baseline, n_recruit))
      phagocyte_activity_engulf = c(phagocyte_activity_engulf, rep(activity_engulf_M0_baseline, n_recruit))
      phagocyte_active_age = c(phagocyte_active_age, rep(0, n_recruit))
      
      new_registry = matrix(0, nrow = n_recruit, ncol = cc_phagocyte)
      phagocyte_bacteria_registry = rbind(phagocyte_bacteria_registry, new_registry)
    }
  }
  
  # Recruit from RIGHT border (x = grid_size, y > 1 to avoid epithelium)
  for (y in 2:grid_size) {
    danger_at_border = danger_signal_grid[y, grid_size]
    n_recruit = rpois(1, lambda = recruitment_rate_danger*danger_at_border)
    
    if (n_recruit > 0) {
      phagocyte_x = c(phagocyte_x, rep(grid_size, n_recruit))
      phagocyte_y = c(phagocyte_y, rep(y, n_recruit))
      phagocyte_pathogens_engulfed = c(phagocyte_pathogens_engulfed, rep(0, n_recruit))
      phagocyte_commensals_engulfed = c(phagocyte_commensals_engulfed, rep(0, n_recruit))
      phagocyte_num_times_activated = c(phagocyte_num_times_activated, rep(0, n_recruit))
      phagocyte_phenotype = c(phagocyte_phenotype, rep(0, n_recruit))
      phagocyte_activity_ROS = c(phagocyte_activity_ROS, rep(activity_ROS_M0_baseline, n_recruit))
      phagocyte_activity_engulf = c(phagocyte_activity_engulf, rep(activity_engulf_M0_baseline, n_recruit))
      phagocyte_active_age = c(phagocyte_active_age, rep(0, n_recruit))
      
      new_registry = matrix(0, nrow = n_recruit, ncol = cc_phagocyte)
      phagocyte_bacteria_registry = rbind(phagocyte_bacteria_registry, new_registry)
    }
  }
}

# ========================================================================
# MOVE PHAGOCYTES AND TREGS
# ========================================================================
density_matrix_tregs = if (randomize_tregs == 1) 0*DAMPs else DAMPs

# Combine DAMPs and PAMPs for macrophage chemotaxis
# This allows macrophages to respond to both tissue damage (DAMPs) and pathogen signals (PAMPs)
density_matrix_phagocytes = DAMPs + PAMPs

all_equal_treg = all(density_matrix_tregs == density_matrix_tregs[1, 1])
all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])

# Move tregs
if (!all_equal_treg) {
  for (i in 1:length(treg_x)) {
    x = treg_x[i]
    y = treg_y[i]
    x_range = max(1, x - 1):min(grid_size, x + 1)
    y_range = max(1, y - 1):min(grid_size, y + 1)
    neighbors_x = rep(x_range, each = length(y_range))
    neighbors_y = rep(y_range, times = length(x_range))
    neighbor_densities = density_matrix_tregs[cbind(neighbors_y, neighbors_x)]
    total = sum(neighbor_densities)
    if (total > 0) {
      probs = neighbor_densities / total
    } else {
      probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
    }
    chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
    treg_x[i] = neighbors_x[chosen_idx]
    treg_y[i] = neighbors_y[chosen_idx]
  }
} else {
  dy_treg = ifelse(treg_y == 1,
                   sample(c(1), size = length(treg_y), replace = TRUE),
                   sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
  dx_treg = iszero_coordinates(dy_treg)
  treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
  treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
}

# Move phagocytes
if (!all_equal_phagocytes) {
  for (i in 1:length(phagocyte_x)) {
    x = phagocyte_x[i]
    y = phagocyte_y[i]
    x_range = max(1, x - 1):min(grid_size, x + 1)
    y_range = max(1, y - 1):min(grid_size, y + 1)
    neighbors_x = rep(x_range, each = length(y_range))
    neighbors_y = rep(y_range, times = length(x_range))
    neighbor_densities = density_matrix_phagocytes[cbind(neighbors_y, neighbors_x)]
    total = sum(neighbor_densities)
    if (total > 0) {
      probs = neighbor_densities / total
    } else {
      probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
    }
    chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
    phagocyte_x[i] = neighbors_x[chosen_idx]
    phagocyte_y[i] = neighbors_y[chosen_idx]
  }
} else {
  dy_phagocyte = ifelse(phagocyte_y == 1,
                        sample(c(1), size = length(phagocyte_y), replace = TRUE),
                        sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
  dx_phagocyte = iszero_coordinates(dy_phagocyte)
  phagocyte_x = pmin(pmax(phagocyte_x + dx_phagocyte, 1), grid_size)
  phagocyte_y = pmin(pmax(phagocyte_y + dy_phagocyte, 1), grid_size)
}

# ========================================================================
# ADD NEW MICROBES
# ========================================================================
n_pathogens_lp_new = round(mean(epithelium$level_injury)*rate_leak_pathogen_injury *
                             length(injury_site_updated))
if (n_pathogens_lp_new > 0) {
  new_pathogen_coords = matrix(c(
    sample(1:grid_size, n_pathogens_lp_new, replace = TRUE, prob = epithelium$level_injury),
    rep(1, n_pathogens_lp_new),
    last_id_pathogen + seq(1, n_pathogens_lp_new)
  ), ncol = 3)
  colnames(new_pathogen_coords) = c("x", "y", "id")
  pathogen_coords = rbind(pathogen_coords, new_pathogen_coords)
  last_id_pathogen = last_id_pathogen + n_pathogens_lp_new
}

n_commensals_lp_new_injury = round(mean(epithelium$level_injury)*rate_leak_commensal_injury *
                                     length(injury_site_updated))
n_commensals_lp_new_baseline = round(rate_leak_commensal_baseline*grid_size)

total_new_commensals = n_commensals_lp_new_baseline + n_commensals_lp_new_injury
if (total_new_commensals > 0) {
  baseline_x = sample(1:grid_size, n_commensals_lp_new_baseline, TRUE)
  injury_x = if (n_commensals_lp_new_injury > 0) {
    sample(1:grid_size, n_commensals_lp_new_injury, TRUE, prob = epithelium$level_injury)
  } else {
    numeric(0)
  }
  new_commensal_coords = matrix(c(
    c(baseline_x, injury_x),
    rep(1, total_new_commensals),
    last_id_commensal + seq(1, total_new_commensals)
  ), ncol = 3)
  colnames(new_commensal_coords) = c("x", "y", "id")
  commensal_coords = rbind(commensal_coords, new_commensal_coords)
  last_id_commensal = last_id_commensal + total_new_commensals
}

# ========================================================================
# PLOTTING
# ========================================================================
if(exists("plot_on")){
  if (plot_on == 1 & (t %% plot_every == 0 | t == 1)) {
    source('./MISC/CONVERT_TO_DATAFRAME_ABM.R')
    p = plot_simtime_simple()
    ggsave(
      paste0(dir_name_frames, "/frame_param_", param_set_id_use, "_rep_", reps_in,
             "_ind_", scenario_ind, "_", t, ".png"),
      plot = p,
      width = 12,
      height = 10,
      dpi = 600,
      bg = "white"
    )
  }
}

# ========================================================================
# UPDATE PHAGOCYTE PHENOTYPES
# C++ ACCELERATION: Batch signal calculations
# ========================================================================
M0_indices = which(phagocyte_phenotype == 0)
M1_indices = which(phagocyte_phenotype == 1)
M2_indices = which(phagocyte_phenotype == 2)

# Process M0 phagocytes (candidates for activation)
if (length(M0_indices) > 0) {
  signals = calculate_phagocyte_signals_cpp(
    M0_indices, phagocyte_x, phagocyte_y, # phagocyte_bacteria_registry,
    act_radius_DAMPs, act_radius_SAMPs, act_radius_PAMPs,
    DAMPs, SAMPs, PAMPs, grid_size
  )
  avg_DAMPs_vec = signals$avg_DAMPs
  avg_SAMPs_vec = signals$avg_SAMPs
  avg_PAMPs_vec = signals$avg_PAMPs  # First return value contains PAMPs average
  # bacteria_count_vec = signals$bacteria_counts

  for (idx in seq_along(M0_indices)) {
    i = M0_indices[idx]
    avg_DAMPs = avg_DAMPs_vec[idx]
    avg_SAMPs = avg_SAMPs_vec[idx]
    avg_PAMPs = avg_PAMPs_vec[idx]
    # bacteria_count = bacteria_count_vec[idx]
    
    # Combine DAMPs + PAMPs as danger signal
    danger_signal      = avg_DAMPs + avg_PAMPs
    
    SAMPS_above_th     = avg_SAMPs >= activation_threshold_SAMPs
    SAMPS_below_th     = avg_SAMPs < activation_threshold_SAMPs
    
    DANGER_above_th    = danger_signal >= activation_threshold_danger
    DANGER_below_th    = danger_signal < activation_threshold_danger
    
    SAMPS_above_DANGER = avg_SAMPs > danger_signal
    DANGER_above_SAMPS = avg_SAMPs <= danger_signal
    
    SAMPS_diff  = max(0,avg_SAMPs-activation_threshold_SAMPs)/activation_threshold_SAMPs
    DANGER_diff = max(0,danger_signal-activation_threshold_danger)/activation_threshold_danger
    
    # why two conditions? because both can be very low, then any activation shouldn't happen
    # and both can be very high, exceeding the thresholds - but then activation should 
    # depend on the one that is dominating.
    if(SAMPS_diff>0 | DANGER_diff>0){
      if (DANGER_diff>SAMPS_diff) {
        # if (DANGER_above_th && DANGER_above_SAMPS) {
        # if (DANGER_above_th) {
        phagocyte_phenotype[i] = 1
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M1_baseline # + activity_ROS_M1_step*bacteria_count
        phagocyte_activity_engulf[i] = activity_engulf_M1_baseline # + activity_engulf_M1_step*bacteria_count
      }else{
        # ## Should I really prevent M0 → M2 transitions so that Tregs can only suppress already-activated M1 macrophages (M1 → M2), not activate naive M0 macrophages.
        # phagocyte_phenotype[i] = 2
        # phagocyte_active_age[i] = 1
        # phagocyte_activity_ROS[i] = activity_ROS_M2_baseline # + activity_ROS_M1_step*bacteria_count
        # phagocyte_activity_engulf[i] = activity_engulf_M2_baseline # + activity_engulf_M1_step*bacteria_count
      }
      ## Prevent M0 → M2 transitions so that Tregs can only suppress already-activated M1 macrophages (M1 → M2), not activate naive M0 macrophages.
      # else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > danger_signal) {
      #   phagocyte_phenotype[i] = 2
      #   phagocyte_active_age[i] = 1
      #   phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
      #   phagocyte_activity_engulf[i] = activity_engulf_M2_baseline # + activity_engulf_M2_step*bacteria_count
      # }
    }
    
  }
}

# Process M1/M2 phagocytes
active_indices = c(M1_indices, M2_indices)
if (length(active_indices) > 0) {
  phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
  old_enough = phagocyte_active_age[active_indices] >= active_age_limit
  candidates = active_indices[old_enough]
  
  if (length(candidates) > 0) {
    # C++ ACCELERATION: Calculate all signals at once
    signals = calculate_phagocyte_signals_cpp(
      candidates, phagocyte_x, phagocyte_y, # phagocyte_bacteria_registry,
      act_radius_DAMPs, act_radius_SAMPs, act_radius_PAMPs,
      DAMPs, SAMPs, PAMPs, grid_size
    )
    avg_DAMPs_vec = signals$avg_DAMPs
    avg_SAMPs_vec = signals$avg_SAMPs
    avg_PAMPs_vec = signals$avg_PAMPs
    # bacteria_count_vec = signals$bacteria_counts # WILL ALWAYS RETURN 0
    
    for (idx in seq_along(candidates)) {
      i = candidates[idx]
      avg_DAMPs = avg_DAMPs_vec[idx]
      avg_SAMPs = avg_SAMPs_vec[idx]
      avg_PAMPs = avg_PAMPs_vec[idx]
      # bacteria_count = bacteria_count_vec[idx]
      
      # Combine DAMPs + PAMPs as danger signal
      danger_signal = avg_DAMPs + avg_PAMPs
      
      # FIX: Phenotype-specific deactivation logic
      # - M1 should deactivate when danger drops (ROS production should stop)
      # - M2 can be maintained by SAMPs (tolerance to commensals)
      # - But if BOTH signals are low, always deactivate
      
      current_phenotype = phagocyte_phenotype[i]
      
      SAMPS_above_th     = avg_SAMPs >= activation_threshold_SAMPs
      SAMPS_below_th     = avg_SAMPs < activation_threshold_SAMPs
      
      DANGER_above_th    = danger_signal >= activation_threshold_danger
      DANGER_below_th    = danger_signal < activation_threshold_danger
      
      SAMPS_above_DANGER = avg_SAMPs > danger_signal
      DANGER_above_SAMPS = avg_SAMPs <= danger_signal
      
      SAMPS_diff  = max(0,avg_SAMPs-activation_threshold_SAMPs)/activation_threshold_SAMPs
      DANGER_diff = max(0,danger_signal-activation_threshold_danger)/activation_threshold_danger
      
      if(current_phenotype == 1){
        # if (SAMPS_below_th && DANGER_below_th) {
        if (SAMPS_diff==0 && DANGER_diff==0) {
          # M1 → M0
          # Both signals low → always deactivate REGARDLESS of phenotype
          phagocyte_phenotype[i] = 0
          phagocyte_active_age[i] = 0 # active age 0 for M0, they don't stay in this phenotype!
          phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          # }else if (SAMPS_above_th && SAMPS_above_DANGER) {
        }else if (SAMPS_diff>DANGER_diff) {
          # M1 → M2
          # SAMPs dominant → M2
          phagocyte_phenotype[i] = 2
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M2_baseline
        }
        #else {
        #   # M1 → M1
        #   # REINFORCE → Reset age given stimulus
        #   phagocyte_phenotype[i] = 1
        #   phagocyte_active_age[i] = 1
        #   phagocyte_activity_ROS[i] = activity_ROS_M1_baseline
        #   phagocyte_activity_engulf[i] = activity_engulf_M1_baseline
        # }
      }else if(current_phenotype == 2){
        # if (SAMPS_below_th && DANGER_below_th) {
        if (SAMPS_diff==0 && DANGER_diff==0) {
          # M2 → M0
          # Both signals low → always deactivate REGARDLESS of phenotype
          phagocyte_phenotype[i] = 0
          phagocyte_active_age[i] = 0
          phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          # }else if (DANGER_above_th && DANGER_above_SAMPS) {
        }else if (SAMPS_diff<DANGER_diff) {
          # M2 → M1
          # DANGER dominant → M1
          phagocyte_phenotype[i] = 1
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M1_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M1_baseline
        }
        # else{
        #   # M2 → M2
        #   # REINFORCE → Reset age given stimulus
        #   phagocyte_phenotype[i] = 2
        #   phagocyte_active_age[i] = 1
        #   phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
        #   phagocyte_activity_engulf[i] = activity_engulf_M2_baseline
        # }
      }
    }
  }
}

# ========================================================================
# UPDATE TREG ACTIVE AGE
# ========================================================================
active_treg_indices = which(treg_phenotype == 1)
if (length(active_treg_indices) > 0) {
  old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
  young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]
  
  if (length(young_tregs) > 0) {
    treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
  }
  
  if (length(old_tregs) > 0) {
    treg_phenotype[old_tregs] = 0
    treg_active_age[old_tregs] = 0
    treg_activity_SAMPs_binary[old_tregs] = 0
  }
}

# ========================================================================
# ENGULFMENT PROCESS
# ========================================================================
phagocyte_positions = paste(phagocyte_x, phagocyte_y, sep = "_")

phagocytes_that_engulfed = rep(FALSE, length(phagocyte_x))

for (i in 1:length(phagocyte_x)) {
  px = phagocyte_x[i]
  py = phagocyte_y[i]
  
  # Pathogen engulfment
  if (nrow(pathogen_coords) > 0) {
    pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
    pathogen_indices = which(pathogen_overlap)
    
    if (length(pathogen_indices) > 0) {
      engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
      indices_to_engulf = pathogen_indices[engulf_success]
      
      if (length(indices_to_engulf) > 0) {
        pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
        
        # C++ ACCELERATION: shift_insert with -1 for pathogens
        phagocyte_bacteria_registry[i, ] = shift_insert_fast_cpp(
          phagocyte_bacteria_registry[i, ],
          rep(-1, length(indices_to_engulf))
        )
        
        phagocytes_that_engulfed[i] = TRUE
        
        phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
        pathogens_killed_by_Mac[phagocyte_phenotype_index] =
          pathogens_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
      }
    }
  }
  
  # Commensal engulfment
  if (nrow(commensal_coords) > 0) {
    commensal_overlap = (commensal_coords[, "x"] == px) & (commensal_coords[, "y"] == py)
    commensal_indices = which(commensal_overlap)
    
    if (length(commensal_indices) > 0) {
      engulf_success = runif(length(commensal_indices)) < phagocyte_activity_engulf[i]
      indices_to_engulf = commensal_indices[engulf_success]
      
      if (length(indices_to_engulf) > 0) {
        commensal_coords = commensal_coords[-indices_to_engulf, , drop = FALSE]
        
        # C++ ACCELERATION: shift_insert with +1 for commensals
        phagocyte_bacteria_registry[i, ] = shift_insert_fast_cpp(
          phagocyte_bacteria_registry[i, ],
          rep(1, length(indices_to_engulf))
        )
        
        phagocytes_that_engulfed[i] = TRUE
        
        phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
        commensals_killed_by_Mac[phagocyte_phenotype_index] =
          commensals_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
      }
    }
  }
  
}

# ========================================================================
# SHIFT REGISTRY FOR PHAGOCYTES THAT DIDN'T ENGULF (AGING MEMORY)
# ========================================================================
phagocytes_to_shift = which(!phagocytes_that_engulfed)
if (length(phagocytes_to_shift) > 0) {
  for (i in phagocytes_to_shift) {
    # Shift registry without adding new bacteria (insert 0 = empty)
    phagocyte_bacteria_registry[i, ] = shift_insert_fast_cpp(
      phagocyte_bacteria_registry[i, ],
      numeric(0)  # Empty vector = shift only, no insertion
    )
  }
}

# ========================================================================
# CALCULATE COUNTS FROM REGISTRY
# ========================================================================
# Count commensals (+1) and pathogens (-1) from registry for each phagocyte
phagocyte_commensals_engulfed = rowSums(phagocyte_bacteria_registry > 0)
phagocyte_pathogens_engulfed  = rowSums(phagocyte_bacteria_registry < 0)

# ========================================================================
# TREG ACTIVATION & EFFECTOR ACTIONS
# C++ ACCELERATION: find_nearby_tregs
# ========================================================================
M1_phagocyte_indices = which(phagocyte_phenotype == 1)
M2_phagocyte_indices = which(phagocyte_phenotype == 2)
M_activate_phagocyte_indices = c(M1_phagocyte_indices, M2_phagocyte_indices)

if (allow_tregs == 1 && length(M_activate_phagocyte_indices) > 0) {
  for (i in M_activate_phagocyte_indices) {
    px = phagocyte_x[i]
    py = phagocyte_y[i]
    
    # C++ ACCELERATION: find nearby tregs
    nearby_treg_indices = find_nearby_tregs_cpp(
      px, py, treg_x, treg_y, treg_vicinity_effect
    )
    
    if (length(nearby_treg_indices) > 0) {
      num_pat_antigens = phagocyte_pathogens_engulfed[i]
      num_com_antigens = phagocyte_commensals_engulfed[i]
      
      if ((num_pat_antigens + num_com_antigens) > 0) {
        rat_com_pat_real = num_com_antigens / (num_com_antigens + num_pat_antigens)# rat_com_pat_real can be interpreted 
        # as the frequency of presenting a commensal antigen to the Treg among all the engulfed antigens
        
        # Step 1: Which antigen type is presented? (based on actual composition)
        commensal_presented = runif(1) < rat_com_pat_real
        
        # Step 2: Does the Treg correctly identify what was presented?
        if (commensal_presented) {
          # A commensal antigen was presented
          # Treg correctly identifies it as commensal with probability = efficiency
          treg_identifies_as_commensal = runif(1) < treg_discrimination_efficiency
        } else {
          # A pathogen antigen was presented
          # Treg incorrectly identifies it as commensal with probability = (1-efficiency)
          treg_identifies_as_commensal = runif(1) < (1 - treg_discrimination_efficiency)
        }
        
        # Step 3: Activate if Treg thinks it saw a commensal
        if (treg_identifies_as_commensal) {
          treg_phenotype[nearby_treg_indices] = 1
          treg_activity_SAMPs_binary[nearby_treg_indices] = 1
          treg_active_age[nearby_treg_indices] = 1
        }
      }
    }
  }
}

# ========================================================================
# KILL MICROBES WITH ROS
# C++ ACCELERATION: Batch killing with single function call (HUGE SPEEDUP)
# ========================================================================
if (nrow(pathogen_coords) > 0) {
  result = kill_microbes_with_ros_cpp(
    pathogen_coords, ROS, act_radius_ROS, th_ROS_microbe, grid_size
  )
  pathogen_coords = result$surviving_microbes
  pathogens_killed_by_ROS = pathogens_killed_by_ROS + result$n_killed
}

if (nrow(commensal_coords) > 0) {
  result = kill_microbes_with_ros_cpp(
    commensal_coords, ROS, act_radius_ROS, th_ROS_microbe, grid_size
  )
  commensal_coords = result$surviving_microbes
  commensals_killed_by_ROS = commensals_killed_by_ROS + result$n_killed
}

# ========================================================================
# UPDATE EPITHELIAL INJURY
# ========================================================================
epithelium_x = epithelium$x

ros_means = calculate_epithelial_ros_cpp(
  epithelium_x, act_radius_ROS, ROS, grid_size
)

# Vectorized injury updates
epithelium$level_injury = epithelium$level_injury +
  logistic_scaled_0_to_5_quantized(pathogen_epithelium_counts)

# ROS-based injury (vectorized)
epithelium$level_injury = epithelium$level_injury + as.integer(ros_means > th_ROS_epith_injury)

# Apply maximum injury constraint (vectorized)
epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)

# RECOVERY: Stochastic recovery (must loop for random number consumption order)
for (i in 1:nrow(epithelium)) {
  if (epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
    epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
  }
}

# ========================================================================
# SAVE ABUNDANCES
# ========================================================================
epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))

phagocyte_counts = c(
  sum(phagocyte_phenotype == 0),
  sum(phagocyte_phenotype == 1),
  sum(phagocyte_phenotype == 2)
)
macrophages_longitudinal[t, ] = phagocyte_counts

microbes_longitudinal[t, ] = c(nrow(commensal_coords), nrow(pathogen_coords))
tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))
microbes_cumdeath_longitudinal[t, ] = c(
  commensals_killed_by_ROS, commensals_killed_by_Mac,
  pathogens_killed_by_ROS, pathogens_killed_by_Mac
)

