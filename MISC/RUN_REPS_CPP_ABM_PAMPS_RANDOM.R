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

evaluation_interval = 300  # Run each random sample for 300 timesteps
max_iterations = 1000      # Number of random samples to try

# Tracking
random_search_results = data.frame(
  iter = integer(),
  objective = numeric(),
  min_epithelial = numeric(),
  mean_epithelial = numeric(),
  min_pathogen = numeric(),
  mean_pathogen = numeric(),
  outcome = character(),  # "success", "collapse", "incomplete"
  termination_time = integer(),
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
cat(sprintf("Evaluation interval: %d timesteps\n", evaluation_interval))
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

  # Run simulation for evaluation_interval timesteps
  score_window_e = numeric(0)
  score_window_p = numeric(0)
  outcome = "incomplete"
  termination_time = evaluation_interval

  for (t in 1:evaluation_interval) {

    # Calculate current scores
    current_epithelial_score = 6*sum(epithelium$level_injury == 0) +
      5*sum(epithelium$level_injury == 1) +
      4*sum(epithelium$level_injury == 2) +
      3*sum(epithelium$level_injury == 3) +
      2*sum(epithelium$level_injury == 4) +
      1*sum(epithelium$level_injury == 5)

    current_pathogens = nrow(pathogen_coords)

    score_window_e = c(score_window_e, current_epithelial_score)
    score_window_p = c(score_window_p, current_pathogens)

    # Check for COLLAPSE
    if (length(score_window_e) >= collapse_duration) {
      recent_scores_collapse = tail(score_window_e, collapse_duration)
      pct_below_threshold = sum(recent_scores_collapse < collapse_threshold) / length(recent_scores_collapse)

      if (pct_below_threshold > collapse_rate) {
        outcome = "collapse"
        termination_time = t
        cat(sprintf("    [COLLAPSE at t=%d] Epithelial score too low\n", t))
        break
      }
    }

    # Check for SUCCESS
    if (length(score_window_e) >= success_duration) {
      recent_scores_success_e = tail(score_window_e, success_duration)
      recent_scores_success_p = tail(score_window_p, success_duration)

      pct_above_threshold_e = sum(recent_scores_success_e > success_threshold_e) / length(recent_scores_success_e)
      pct_above_threshold_p = sum(recent_scores_success_p < success_threshold_p) / length(recent_scores_success_p)

      is_success_e = (pct_above_threshold_e > success_rate)
      is_success_p = (pct_above_threshold_p > success_rate)

      if (is_success_e && is_success_p) {
        outcome = "success"
        termination_time = t
        cat(sprintf("    [SUCCESS at t=%d] Criteria met!\n", t))

        # Log success
        success_line = sprintf("iter=%d | t=%d | score=%.1f | theta=[%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                               iter, t, min(score_window_e),
                               current_theta[1], current_theta[2], current_theta[3],
                               current_theta[4], current_theta[5])
        cat(success_line, file = success_log_file, append = TRUE)
        break
      }
    }

    # Run one simulation step
    source("./MISC/MAIN_SIM.R")
  }

  # Calculate objective
  objective = min(score_window_e)

  cat(sprintf("    Outcome: %s | Objective: %.1f (min_epithelial) | Time: %d timesteps\n",
              outcome, objective, termination_time))
  cat(sprintf("    Epithelial - min: %.1f, mean: %.1f\n",
              min(score_window_e), mean(score_window_e)))
  cat(sprintf("    Pathogens - min: %d, mean: %.1f\n\n",
              min(score_window_p), mean(score_window_p)))

  # Record results
  random_search_results = rbind(random_search_results, data.frame(
    iter = iter,
    objective = objective,
    min_epithelial = min(score_window_e),
    mean_epithelial = mean(score_window_e),
    min_pathogen = min(score_window_p),
    mean_pathogen = mean(score_window_p),
    outcome = outcome,
    termination_time = termination_time,
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
cat(sprintf("Successes: %d\n", sum(random_search_results$outcome == "success")))
cat(sprintf("Collapses: %d\n", sum(random_search_results$outcome == "collapse")))
cat(sprintf("Incomplete: %d\n", sum(random_search_results$outcome == "incomplete")))

if (sum(random_search_results$outcome == "success") > 0) {
  best_success = random_search_results[random_search_results$outcome == "success", ]
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
  cat(sprintf("  Iter: %d | Objective: %.1f | Outcome: %s\n",
              overall_best$iter, overall_best$objective, overall_best$outcome))
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
