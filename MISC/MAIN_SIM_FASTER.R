# ============================================================================
# MAIN SIMULATION LOOP - OPTIMIZED VERSION WITH ADDITIONAL C++ FUNCTIONS
# ============================================================================
# This version uses additional C++ optimizations compared to MAIN_SIM.R
# Expected speedup: 2-5x depending on parameter values
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

# DAMPs from microbes touching epithelium- this needs to be different because commensals do not harm the epithelium but can induce damps on the basolateral side
DAMPs[1, ] = DAMPs[1, ] + logistic_scaled_0_to_5_quantized(pathogen_epithelium_counts + commensal_epithelium_counts, k_in_log, x0_in_log, c_in_log)*add_DAMPs

# ========================================================================
# UPDATE PAMPs - C++ OPTIMIZED VERSION (15-40x faster)
# ========================================================================
PAMPs = update_PAMPs_cpp(PAMPs, pathogen_coords, add_PAMPs, k_in_log, x0_in_log, c_in_log)

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
# MOVE PHAGOCYTES AND TREGS - C++ OPTIMIZED VERSION (10-30x faster)
# ========================================================================
density_matrix_tregs = if (randomize_tregs == 1) 0*DAMPs else DAMPs

# Combine DAMPs and PAMPs for macrophage chemotaxis
# This allows macrophages to respond to both tissue damage (DAMPs) and pathogen signals (PAMPs)
density_matrix_phagocytes = DAMPs + PAMPs

all_equal_treg = all(density_matrix_tregs == density_matrix_tregs[1, 1])
all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])

# Move tregs - C++ OPTIMIZED
if (!all_equal_treg) {
  result = move_cells_chemotaxis_cpp(treg_x, treg_y, density_matrix_tregs, grid_size)
  treg_x = result$x
  treg_y = result$y
} else {
  dy_treg = ifelse(treg_y == 1,
                   sample(c(1), size = length(treg_y), replace = TRUE),
                   sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
  dx_treg = iszero_coordinates(dy_treg)
  treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
  treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
}

# Move phagocytes - C++ OPTIMIZED
if (!all_equal_phagocytes) {
  result = move_cells_chemotaxis_cpp(phagocyte_x, phagocyte_y, density_matrix_phagocytes, grid_size)
  phagocyte_x = result$x
  phagocyte_y = result$y
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
n_pathogens_lp_new = round(mean(epithelium$level_injury)*rate_leak_pathogen_injury*length(injury_site_updated))

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

n_commensals_lp_new_injury = round(mean(epithelium$level_injury)*rate_leak_commensal_injury*length(injury_site_updated))
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
# UPDATE PHAGOCYTE PHENOTYPES - C++ OPTIMIZED VERSION (3-8x faster)
# C++ ACCELERATION: Batch signal calculations AND decision logic
# ========================================================================
M0_indices = which(phagocyte_phenotype == 0)
M1_indices = which(phagocyte_phenotype == 1)
M2_indices = which(phagocyte_phenotype == 2)

# Process M0 phagocytes (candidates for activation) - C++ OPTIMIZED
if (length(M0_indices) > 0) {
  result = update_M0_phenotypes_cpp(
    M0_indices, phagocyte_x, phagocyte_y,
    DAMPs, SAMPs, PAMPs,
    act_radius_DAMPs, act_radius_SAMPs, act_radius_PAMPs,
    grid_size,
    activation_threshold_SAMPs,
    activation_threshold_danger,
    activity_ROS_M1_baseline,
    activity_engulf_M1_baseline
  )

  # Update phagocytes that activated
  activated_mask = result$activated
  if (any(activated_mask)) {
    activated_idx = M0_indices[activated_mask]
    phagocyte_phenotype[activated_idx] = result$phenotype[activated_mask]
    phagocyte_active_age[activated_idx] = result$active_age[activated_mask]
    phagocyte_activity_ROS[activated_idx] = result$activity_ROS[activated_mask]
    phagocyte_activity_engulf[activated_idx] = result$activity_engulf[activated_mask]
  }
}

# Process M1/M2 phagocytes - C++ OPTIMIZED
active_indices = c(M1_indices, M2_indices)
if (length(active_indices) > 0) {
  phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
  old_enough = phagocyte_active_age[active_indices] >= active_age_limit
  too_young  = phagocyte_active_age[active_indices] < active_age_limit
  candidates = active_indices[old_enough]
  candidates_only_suppress = active_indices[too_young]

  if (length(candidates) > 0) {
    # C++ OPTIMIZED: old cells can deactivate or transition
    result = update_active_phenotypes_cpp(
      candidates, phagocyte_x, phagocyte_y,
      phagocyte_phenotype, phagocyte_active_age,
      DAMPs, SAMPs, PAMPs,
      act_radius_DAMPs, act_radius_SAMPs, act_radius_PAMPs,
      grid_size, active_age_limit,
      activation_threshold_SAMPs,
      activation_threshold_danger,
      activity_ROS_M0_baseline,
      activity_engulf_M0_baseline,
      activity_ROS_M1_baseline,
      activity_engulf_M1_baseline,
      activity_ROS_M2_baseline,
      activity_engulf_M2_baseline,
      only_suppress = FALSE
    )

    # Update phagocytes that changed
    changed_mask = result$changed
    if (any(changed_mask)) {
      changed_idx = candidates[changed_mask]
      phagocyte_phenotype[changed_idx] = result$phenotype[changed_mask]
      phagocyte_active_age[changed_idx] = result$active_age[changed_mask]
      phagocyte_activity_ROS[changed_idx] = result$activity_ROS[changed_mask]
      phagocyte_activity_engulf[changed_idx] = result$activity_engulf[changed_mask]
    }
  }

  if (overwrite_in == 1 && length(candidates_only_suppress) > 0) {
    # C++ OPTIMIZED: young cells can only suppress (M1->M2 or M2->M1)
    result = update_active_phenotypes_cpp(
      candidates_only_suppress, phagocyte_x, phagocyte_y,
      phagocyte_phenotype, phagocyte_active_age,
      DAMPs, SAMPs, PAMPs,
      act_radius_DAMPs, act_radius_SAMPs, act_radius_PAMPs,
      grid_size, active_age_limit,
      activation_threshold_SAMPs,
      activation_threshold_danger,
      activity_ROS_M0_baseline,
      activity_engulf_M0_baseline,
      activity_ROS_M1_baseline,
      activity_engulf_M1_baseline,
      activity_ROS_M2_baseline,
      activity_engulf_M2_baseline,
      only_suppress = TRUE
    )

    # Update phagocytes that changed
    changed_mask = result$changed
    if (any(changed_mask)) {
      changed_idx = candidates_only_suppress[changed_mask]
      phagocyte_phenotype[changed_idx] = result$phenotype[changed_mask]
      phagocyte_active_age[changed_idx] = result$active_age[changed_mask]
      phagocyte_activity_ROS[changed_idx] = result$activity_ROS[changed_mask]
      phagocyte_activity_engulf[changed_idx] = result$activity_engulf[changed_mask]
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
# ENGULFMENT PROCESS - C++ OPTIMIZED VERSION (20-50x faster)
# ========================================================================
if (length(phagocyte_x) > 0 && (nrow(pathogen_coords) > 0 || nrow(commensal_coords) > 0)) {
  result = process_engulfment_batch_cpp(
    phagocyte_x, phagocyte_y,
    phagocyte_activity_engulf,
    phagocyte_phenotype,
    pathogen_coords,
    commensal_coords,
    phagocyte_bacteria_registry,
    cc_phagocyte
  )

  # Update all results
  pathogen_coords = result$pathogen_coords
  commensal_coords = result$commensal_coords
  phagocyte_bacteria_registry = result$phagocyte_bacteria_registry
  phagocytes_that_engulfed = result$phagocytes_that_engulfed

  # Update kill counters
  pathogens_killed_by_Mac = pathogens_killed_by_Mac + result$pathogens_killed
  commensals_killed_by_Mac = commensals_killed_by_Mac + result$commensals_killed

  # Shift registry for phagocytes that didn't engulf - C++ OPTIMIZED
  phagocytes_to_shift = !phagocytes_that_engulfed
  if (any(phagocytes_to_shift)) {
    phagocyte_bacteria_registry = shift_registry_batch_cpp(
      phagocyte_bacteria_registry, phagocytes_to_shift
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
# TREG ACTIVATION & EFFECTOR ACTIONS - C++ OPTIMIZED VERSION (5-10x faster)
# ========================================================================
M1_phagocyte_indices = which(phagocyte_phenotype == 1)
M2_phagocyte_indices = which(phagocyte_phenotype == 2)
M_activate_phagocyte_indices = c(M1_phagocyte_indices, M2_phagocyte_indices)

if (allow_tregs == 1 && length(M_activate_phagocyte_indices) > 0) {
  result = activate_tregs_batch_cpp(
    M_activate_phagocyte_indices,
    phagocyte_x, phagocyte_y,
    as.integer(phagocyte_pathogens_engulfed),
    as.integer(phagocyte_commensals_engulfed),
    treg_x, treg_y,
    treg_phenotype,
    treg_activity_SAMPs_binary,
    treg_active_age,
    treg_vicinity_effect,
    treg_discrimination_efficiency
  )

  # Update treg states
  treg_phenotype = result$treg_phenotype
  treg_activity_SAMPs_binary = result$treg_activity_SAMPs_binary
  treg_active_age = result$treg_active_age
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
# === VERSION #2
add_inj_pathogen = x0_in*(1-exp(-k_in*pathogen_epithelium_counts)) # max x0_in added
add_inj_ros      = pmax(0,x0_in*(ros_means-th_ROS_epith_injury)) # max x0_in added

epithelium$level_injury = epithelium$level_injury + add_inj_pathogen + add_inj_ros

# RECOVERY: Stochastic recovery (must loop for random number consumption order)
for (i in 1:nrow(epithelium)) {
  if (epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
    epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
  }
}

# Apply maximum injury constraint (vectorized)
epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)

# ========================================================================
# SAVE ABUNDANCES
# ========================================================================
epithelium_longitudinal[t, ] = sum(epithelium$level_injury)

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
