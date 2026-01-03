# ============================================================================
# C++ ACCELERATED VERSION WITH PAMPs SUPPORT + RANDOM SEARCH OPTIMIZATION
# ============================================================================

# Helper function to reset simulation state
reset_simulation_state = function() {
  DAMPs <<- matrix(0, grid_size, grid_size)
  SAMPs <<- matrix(0, grid_size, grid_size)
  PAMPs <<- matrix(0, grid_size, grid_size)
  ROS   <<- matrix(0, grid_size, grid_size)

  epithelium$level_injury <<- 0
  epithelium[injury_site, ]$level_injury <<- 1

  if (n_pathogens_lp == 0) {
    pathogen_coords <<- matrix(numeric(0), ncol = 3)
    colnames(pathogen_coords) <<- c("x", "y", "id")
  } else {
    pathogen_coords <<- matrix(c(
      sample(injury_site, n_pathogens_lp, TRUE),
      rep(1, n_pathogens_lp),
      seq(1, n_pathogens_lp)
    ), ncol = 3)
    colnames(pathogen_coords) <<- c("x", "y", "id")
  }

  commensal_coords <<- matrix(c(
    sample(1:grid_size, n_commensals_lp, TRUE),
    sample(1:grid_size, n_commensals_lp, TRUE),
    seq(1, n_commensals_lp)
  ), ncol = 3)
  colnames(commensal_coords) <<- c("x", "y", "id")

  last_id_pathogen <<- n_pathogens_lp
  last_id_commensal <<- n_commensals_lp

  phagocyte_x <<- sample(1:grid_size, n_phagocytes, TRUE)
  phagocyte_y <<- sample(1:grid_size, n_phagocytes, TRUE)
  phagocyte_phenotype <<- rep(0, n_phagocytes)
  phagocyte_active_age <<- rep(0, n_phagocytes)
  phagocyte_activity_ROS <<- rep(activity_ROS_M0_baseline, n_phagocytes)
  phagocyte_activity_engulf <<- rep(activity_engulf_M0_baseline, n_phagocytes)
  phagocyte_bacteria_registry <<- matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)
  phagocyte_pathogens_engulfed <<- rep(0, n_phagocytes)
  phagocyte_commensals_engulfed <<- rep(0, n_phagocytes)
  phagocyte_num_times_activated <<- rep(0, n_phagocytes)

  treg_x <<- sample(1:grid_size, n_tregs, TRUE)
  treg_y <<- sample(1:grid_size, n_tregs, TRUE)
  treg_phenotype <<- rep(0, n_tregs)
  treg_active_age <<- rep(0, n_tregs)
  treg_activity_SAMPs_binary <<- rep(0, n_tregs)
}

# Initialize fields
DAMPs = matrix(0, grid_size, grid_size)
SAMPs = matrix(0, grid_size, grid_size)
PAMPs = matrix(0, grid_size, grid_size)
ROS   = matrix(0, grid_size, grid_size)

# Initialize epithelium
epithelium = data.frame(
  x = seq(1, grid_size, 1),
  y = rep(0, grid_size),
  level_injury = 0,
  id = seq(1, grid_size)
)
epithelium[injury_site, ]$level_injury = 1

# Initialize phagocytes
phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
phagocyte_y = sample(1:grid_size, n_phagocytes, TRUE)
phagocyte_pathogens_engulfed = rep(0, n_phagocytes)
phagocyte_commensals_engulfed = rep(0, n_phagocytes)
phagocyte_num_times_activated = rep(0, n_phagocytes)
phagocyte_phenotype = rep(0, n_phagocytes)
phagocyte_activity_ROS = rep(activity_ROS_M0_baseline, n_phagocytes)
phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
phagocyte_active_age = rep(0, n_phagocytes)
phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)

# Initialize tregs
treg_x = sample(1:grid_size, n_tregs, TRUE)
treg_y = sample(1:grid_size, n_tregs, TRUE)
treg_active_age = rep(0, n_tregs)
treg_phenotype = rep(0, n_tregs)
treg_activity_SAMPs_binary = rep(0, n_tregs)

# Initialize pathogens
if (n_pathogens_lp == 0) {
  pathogen_coords = matrix(numeric(0), ncol = 3)
  colnames(pathogen_coords) = c("x", "y", "id")
} else {
  pathogen_coords = matrix(c(
    sample(injury_site, n_pathogens_lp, TRUE),
    rep(1, n_pathogens_lp),
    seq(1, n_pathogens_lp)
  ), ncol = 3)
  colnames(pathogen_coords) = c("x", "y", "id")
}

# Initialize commensals
commensal_coords = matrix(c(
  sample(1:grid_size, n_commensals_lp, TRUE),
  sample(1:grid_size, n_commensals_lp, TRUE),
  seq(1, n_commensals_lp)
), ncol = 3)
colnames(commensal_coords) = c("x", "y", "id")

last_id_pathogen = n_pathogens_lp
last_id_commensal = n_commensals_lp

# Death counters
pathogens_killed_by_ROS = 0
pathogens_killed_by_Mac = rep(0, 3)
commensals_killed_by_ROS = 0
commensals_killed_by_Mac = rep(0, 3)

# Longitudinal tracking
epithelium_longitudinal  = matrix(0, nrow = t_max, ncol = (max_level_injury + 1))
macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3)
microbes_longitudinal    = matrix(0, nrow = t_max, ncol = 2)
tregs_longitudinal       = matrix(0, nrow = t_max, ncol = 2)
microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2*4)

# ============================================================================
# RANDOM SEARCH CONFIGURATION
# ============================================================================
param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

lower_bounds = param_bounds$lower
upper_bounds = param_bounds$upper

max_iterations = 1000      # Number of random samples to try

# Tracking
random_search_results = data.frame(
  iter = integer(),
  objective = numeric(),
  min_epithelial = numeric(),
  mean_epithelial = numeric(),
  min_pathogen = numeric(),
  mean_pathogen = numeric(),
  is_success = logical(),
  diffusion_speed_SAMPs = numeric(),
  add_SAMPs = numeric(),
  SAMPs_decay = numeric(),
  treg_discrimination_efficiency = numeric(),
  activation_threshold_SAMPs = numeric(),
  stringsAsFactors = FALSE
)

# Initialize success log file
success_log_file = paste0("./random_search_successes_", param_set_id_use, "_scenario_", scenario_ind, "_n2_", n2, ".txt")
cat(paste0("# Random Search Success Log pat_", pat_level, "_ros_", ros_level, "\n"),
    file = success_log_file, append = FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("RANDOM SEARCH OPTIMIZATION STARTING\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Simulation time: %d timesteps\n", t_max))
cat(sprintf("Evaluation window: last %d timesteps\n", success_duration))
cat(sprintf("Max iterations: %d\n", max_iterations))
cat(sprintf("Parameter bounds:\n"))
for (i in 1:length(param_names)) {
  cat(sprintf("  %s: [%.4f, %.4f]\n", param_names[i], lower_bounds[i], upper_bounds[i]))
}
cat(sprintf("========================================\n\n"))

# ============================================================================
# MAIN RANDOM SEARCH LOOP
# ============================================================================
for (iter in 1:max_iterations) {

  # Sample random parameters
  current_theta = sapply(1:length(param_names), function(i) {
    runif(1, lower_bounds[i], upper_bounds[i])
  })

  # Set parameter values
  diffusion_speed_SAMPs          = current_theta[1]
  add_SAMPs                      = current_theta[2]
  SAMPs_decay                    = current_theta[3]
  treg_discrimination_efficiency = current_theta[4]
  activation_threshold_SAMPs     = current_theta[5]

  cat(sprintf("\n>>> [ITERATION %d/%d] <<<\n", iter, max_iterations))
  cat(sprintf("    Testing theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              current_theta[1], current_theta[2], current_theta[3],
              current_theta[4], current_theta[5]))

  # Reset simulation
  reset_simulation_state()

  # Run full simulation (t_max timesteps)
  for (t in 1:t_max) {
    source("./MISC/MAIN_SIM.R")
  }

  # Evaluate ONLY the last success_duration timesteps
  eval_start = max(1, t_max - success_duration + 1)

  eval_scores_e = numeric(success_duration)
  eval_scores_p = numeric(success_duration)

  eval_idx = 1
  for (t_eval in eval_start:t_max) {
    eval_scores_e[eval_idx] = epithelium_longitudinal[t_eval, 1]*6 +
                               epithelium_longitudinal[t_eval, 2]*5 +
                               epithelium_longitudinal[t_eval, 3]*4 +
                               epithelium_longitudinal[t_eval, 4]*3 +
                               epithelium_longitudinal[t_eval, 5]*2 +
                               epithelium_longitudinal[t_eval, 6]*1
    eval_scores_p[eval_idx] = microbes_longitudinal[t_eval, 2]  # pathogen column
    eval_idx = eval_idx + 1
  }

  # Check success criteria on evaluation window
  pct_above_threshold_e = sum(eval_scores_e > success_threshold_e) / length(eval_scores_e)
  pct_above_threshold_p = sum(eval_scores_p < success_threshold_p) / length(eval_scores_p)

  is_success_e = (pct_above_threshold_e > success_rate)
  is_success_p = (pct_above_threshold_p > success_rate)
  is_success = is_success_e & is_success_p

  # Calculate objective (min of evaluation window)
  objective = min(eval_scores_e)

  if (is_success) {
    cat(sprintf("    [SUCCESS] Criteria met in final %d timesteps!\n", success_duration))

    # Log success
    success_line = sprintf("iter=%d | score=%.1f | theta=[%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                           iter, objective,
                           current_theta[1], current_theta[2], current_theta[3],
                           current_theta[4], current_theta[5])
    cat(success_line, file = success_log_file, append = TRUE)
  }

  cat(sprintf("    Success: %s | Objective: %.1f (min_epithelial in last %d steps)\n",
              is_success, objective, success_duration))
  cat(sprintf("    Evaluation window - Epithelial (min: %.1f, mean: %.1f) | Pathogens (min: %d, mean: %.1f)\n\n",
              min(eval_scores_e), mean(eval_scores_e),
              min(eval_scores_p), mean(eval_scores_p)))

  # Record results
  random_search_results = rbind(random_search_results, data.frame(
    iter = iter,
    objective = objective,
    min_epithelial = min(eval_scores_e),
    mean_epithelial = mean(eval_scores_e),
    min_pathogen = min(eval_scores_p),
    mean_pathogen = mean(eval_scores_p),
    is_success = is_success,
    diffusion_speed_SAMPs = current_theta[1],
    add_SAMPs = current_theta[2],
    SAMPs_decay = current_theta[3],
    treg_discrimination_efficiency = current_theta[4],
    activation_threshold_SAMPs = current_theta[5],
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# SAVE RESULTS
# ============================================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("RANDOM SEARCH COMPLETE\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Total iterations: %d\n", max_iterations))
cat(sprintf("Successes: %d\n", sum(random_search_results$is_success)))

if (sum(random_search_results$is_success) > 0) {
  best_success = random_search_results[random_search_results$is_success, ]
  best_success = best_success[which.max(best_success$objective), ]
  cat(sprintf("\nBest successful configuration:\n"))
  cat(sprintf("  Iter: %d | Objective: %.1f\n", best_success$iter, best_success$objective))
  cat(sprintf("  Theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              best_success$diffusion_speed_SAMPs,
              best_success$add_SAMPs,
              best_success$SAMPs_decay,
              best_success$treg_discrimination_efficiency,
              best_success$activation_threshold_SAMPs))
}

if (nrow(random_search_results) > 0) {
  overall_best = random_search_results[which.max(random_search_results$objective), ]
  cat(sprintf("\nOverall best configuration (any outcome):\n"))
  cat(sprintf("  Iter: %d | Objective: %.1f | Success: %s\n",
              overall_best$iter, overall_best$objective, overall_best$is_success))
  cat(sprintf("  Theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              overall_best$diffusion_speed_SAMPs,
              overall_best$add_SAMPs,
              overall_best$SAMPs_decay,
              overall_best$treg_discrimination_efficiency,
              overall_best$activation_threshold_SAMPs))
}

cat(sprintf("========================================\n\n"))

# Save results to CSV
write.csv(random_search_results,
          file = paste0("./random_search_results_", param_set_id_use, "_scenario_", scenario_ind, "_n2_", n2, ".csv"),
          row.names = FALSE)

# Create final longitudinal data (from last simulation)
longitudinal_df = data.frame(
  epithelium_longitudinal,
  macrophages_longitudinal,
  microbes_longitudinal,
  tregs_longitudinal,
  microbes_cumdeath_longitudinal
)

colnames(longitudinal_df) = colnames_insert

longitudinal_df$t = 1:t_max
longitudinal_df$sterile = sterile
longitudinal_df$macspec_on = macspec_on
longitudinal_df$tregs_on = allow_tregs
longitudinal_df$ros_level = ros_level
longitudinal_df$pat_level = pat_level
longitudinal_df$randomize_tregs = randomize_tregs
longitudinal_df$param_set_id = param_set_use$param_set_id

longitudinal_df = longitudinal_df %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy +
                                                      5*epithelial_inj_1 +
                                                      4*epithelial_inj_2 +
                                                      3*epithelial_inj_3 +
                                                      2*epithelial_inj_4 +
                                                      1*epithelial_inj_5)

longitudinal_df$time_ss_e = steady_state_idx(longitudinal_df$epithelial_score)
longitudinal_df$time_ss_p = steady_state_idx(longitudinal_df$pathogen)

longitudinal_df = longitudinal_df %>%
  dplyr::select(t,
                sterile,
                tregs_on,
                ros_level,
                pat_level,
                macspec_on,
                randomize_tregs,
                param_set_id,
                epithelial_score,
                time_ss_e,
                time_ss_p,
                everything())

longitudinal_df_keep = rbind(longitudinal_df_keep, longitudinal_df)
