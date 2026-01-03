# ============================================================================
# C++ ACCELERATED VERSION WITH PAMPs SUPPORT + EMBEDDED SPSA OPTIMIZATION
# ============================================================================
# ============================================================================
# SPSA CONFIGURATION
# ============================================================================
spsa_params = list(
  param_names = c("diffusion_speed_SAMPs", 
                  "add_SAMPs", 
                  "SAMPs_decay", 
                  "treg_discrimination_efficiency",
                  "activation_threshold_SAMPs"),
  
  theta = c(diffusion_speed_SAMPs, 
            add_SAMPs, 
            SAMPs_decay, 
            treg_discrimination_efficiency, 
            activation_threshold_SAMPs),
  
  lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
  upper = c(0.120, 0.500, 0.250, 1.000, 1.000),
  
  c_scale = c(0.005, 0.020, 0.020, 0.05, 0.05),
  a_scale = c(0.001, 0.005, 0.005, 0.01, 0.01),
  
  a = 1.0,
  c = 1.0,
  A = 100,
  alpha = 0.602,
  gamma = 0.101,
  
  k = 0,
  score_history = numeric(0),
  theta_history = matrix(nrow = 0, ncol = 5),
  f_plus_history = numeric(0),
  f_minus_history = numeric(0),
  t_history = numeric(0)
)

# Helper: clip to bounds
clip_theta = function(theta, lower, upper) {
  pmax(pmin(theta, upper), lower)
}

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
# SPSA STATE
# ============================================================================
spsa_update_interval = 100*active_age_limit
spsa_score_window_e = numeric(0)
spsa_score_window_p = numeric(0)
spsa_last_window_mean = NA
spsa_phase = "baseline"
spsa_delta = NULL
spsa_f_plus = NA
spsa_f_minus = NA

# Initialize success log file
cat(paste0("# SPSA Success Log pat_",pat_level,"_ros_",ros_level,"\n"), file = paste0("./spsa_successes_", param_set_id_use,"_scenario_",scenario_ind,"_n2_",n2,".txt"), append = FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("SPSA OPTIMIZATION STARTING\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Update interval: %d timesteps\n", spsa_update_interval))
cat(sprintf("Initial theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
            spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
            spsa_params$theta[4], spsa_params$theta[5]))
cat(sprintf("========================================\n\n"))

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
for (t in 1:t_max) {
  
  # ========================================================================
  # ACCUMULATE EPITHELIAL SCORE FOR SPSA
  # ========================================================================
  current_epithelial_score = 6*sum(epithelium$level_injury == 0) +
    5*sum(epithelium$level_injury == 1) +
    4*sum(epithelium$level_injury == 2) +
    3*sum(epithelium$level_injury == 3) +
    2*sum(epithelium$level_injury == 4) +
    1*sum(epithelium$level_injury == 5)
  
  current_pathogens = nrow(pathogen_coords)
  
  spsa_score_window_e = c(spsa_score_window_e, current_epithelial_score)
  spsa_score_window_p = c(spsa_score_window_p, current_pathogens)
  
  # ========================================================================
  # EARLY DETECTION: COLLAPSE OR SUCCESS (check BEFORE regular update)
  # ========================================================================
  early_triggered = FALSE
  
  # # DEBUG: Print window stats at regular interval boundaries
  # if (spsa_phase != "baseline" && length(spsa_score_window) >= 50 && (t %% spsa_update_interval) == 0) {
  #   cat(sprintf("  DEBUG t=%d: len=%d, mean=%.1f, min=%d, max=%d, all<30=%s, phase=%s\n",
  #               t, length(spsa_score_window), mean(spsa_score_window), 
  #               min(spsa_score_window), max(spsa_score_window),
  #               all(tail(spsa_score_window, 30) < 30), spsa_phase))
  # }
  
  if (spsa_phase != "baseline") {
    
    # Check for COLLAPSE (needs collapse_duration scores)
    is_collapse = FALSE
    if (length(spsa_score_window_e) >= collapse_duration) {
      recent_scores_collapse = tail(spsa_score_window_e, collapse_duration)
      pct_below_threshold = sum(recent_scores_collapse < collapse_threshold) / length(recent_scores_collapse)
      is_collapse = all(recent_scores_collapse < collapse_threshold) ||
        (pct_below_threshold > collapse_rate)
    }
    
    # Check for SUCCESS (needs success_duration scores)
    is_success = FALSE
    if (length(spsa_score_window_e) >= success_duration) {
      recent_scores_success_e = tail(spsa_score_window_e, success_duration)
      recent_scores_success_p = tail(spsa_score_window_p, success_duration)
      
      pct_above_threshold_e = sum(recent_scores_success_e > success_threshold_e) / length(recent_scores_success_e)
      is_success_e = all(recent_scores_success_e > success_threshold_e) || (pct_above_threshold_e > success_rate)
      
      pct_above_threshold_p = sum(recent_scores_success_p < success_threshold_p) / length(recent_scores_success_p)
      is_success_p = all(recent_scores_success_p < success_threshold_p) || (pct_above_threshold_p > success_rate)
      
      is_success = is_success_e & is_success_p
    }
    
    if (is_collapse || is_success) {

      early_triggered = TRUE

      if (is_collapse) {
        early_objective = min(spsa_score_window_e) + min(spsa_score_window_p)
        cat(sprintf("\n>>> [EARLY TERMINATION: COLLAPSE] <<<\n"))
        cat(sprintf("    Time: t=%d | Phase: %s | SPSA iter: %d\n", t, spsa_phase, spsa_params$k))
        cat(sprintf("    Objective: %.1f (min_epithelial + min_pathogen)\n", early_objective))
        cat(sprintf("    Pathogens: %d\n", current_pathogens))
        cat(sprintf("    Old theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))

        # Randomize theta upon collapse to escape bad parameter region
        spsa_params$theta = sapply(1:length(spsa_params$theta), function(i) {
          runif(1, spsa_params$lower[i], spsa_params$upper[i])
        })

        cat(sprintf("    → Randomizing theta to escape bad parameter region\n"))
        cat(sprintf("    New random theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))

        # Update parameter values
        diffusion_speed_SAMPs          = spsa_params$theta[1]
        add_SAMPs                      = spsa_params$theta[2]
        SAMPs_decay                    = spsa_params$theta[3]
        treg_discrimination_efficiency = spsa_params$theta[4]
        activation_threshold_SAMPs     = spsa_params$theta[5]

        # Reset simulation state and restart from baseline
        cat(sprintf("    → Resetting simulation state and restarting from baseline phase\n\n"))
        reset_simulation_state()
        spsa_score_window_e = numeric(0)
        spsa_score_window_p = numeric(0)
        spsa_phase = "baseline"

        next  # Skip the rest of the loop and continue with baseline

      } else {
        early_objective = min(spsa_score_window_e) + min(spsa_score_window_p)
        cat(sprintf("\n>>> [EARLY TERMINATION: SUCCESS] <<<\n"))
        cat(sprintf("    Time: t=%d | Phase: %s | SPSA iter: %d\n", t, spsa_phase, spsa_params$k))
        cat(sprintf("    Objective: %.1f (min_epithelial + min_pathogen)\n", early_objective))
        cat(sprintf("    Pathogens: %d\n", current_pathogens))
        cat(sprintf("    Current theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))

        # Write success to file
        success_line = sprintf("t=%d | iter=%d | score=%.1f | phase=%s | theta=[%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                               t, spsa_params$k, early_objective, spsa_phase,
                               spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                               spsa_params$theta[4], spsa_params$theta[5])

        cat(success_line, file = paste0("./spsa_successes_", param_set_id_use,"_scenario_",scenario_ind,"_n2_",n2,".txt"), append = TRUE)
      }

      if (spsa_phase == "plus") {
        spsa_f_plus = early_objective
        cat(sprintf("    Recorded: f_plus = %.2f\n", spsa_f_plus))
        
        cat(sprintf("    → Resetting simulation state for MINUS evaluation\n"))
        reset_simulation_state()
        spsa_score_window_e = numeric(0)
        spsa_score_window_p = numeric(0)
        
        c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
        perturb = c_k*spsa_params$c_scale*spsa_delta
        
        theta_minus = clip_theta(
          spsa_params$theta - perturb,
          spsa_params$lower, spsa_params$upper
        )
        
        diffusion_speed_SAMPs          = theta_minus[1]
        add_SAMPs                      = theta_minus[2]
        SAMPs_decay                    = theta_minus[3]
        treg_discrimination_efficiency = theta_minus[4]
        activation_threshold_SAMPs     = theta_minus[5]
        
        cat(sprintf("    → Switching to MINUS phase with theta_minus: [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n",
                    theta_minus[1], theta_minus[2], theta_minus[3], theta_minus[4], theta_minus[5]))
        
        spsa_phase = "minus"
        
      } else if (spsa_phase == "minus") {
        spsa_f_minus = early_objective
        cat(sprintf("    Recorded: f_minus = %.2f\n", spsa_f_minus))
        
        c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
        a_k = spsa_params$a / (spsa_params$A + spsa_params$k)^spsa_params$alpha
        
        perturb = c_k*spsa_params$c_scale*spsa_delta
        g_hat = (spsa_f_plus - spsa_f_minus) / (2*perturb)
        
        step = a_k*spsa_params$a_scale*g_hat
        
        spsa_params$theta = clip_theta(
          spsa_params$theta + step,
          spsa_params$lower, spsa_params$upper
        )
        
        diffusion_speed_SAMPs          = spsa_params$theta[1]
        add_SAMPs                      = spsa_params$theta[2]
        SAMPs_decay                    = spsa_params$theta[3]
        treg_discrimination_efficiency = spsa_params$theta[4]
        activation_threshold_SAMPs     = spsa_params$theta[5]
        
        cat(sprintf("\n--- SPSA UPDATE (EARLY, iter %d) ---\n", spsa_params$k))
        cat(sprintf("    f_plus:  %.2f\n", spsa_f_plus))
        cat(sprintf("    f_minus: %.2f\n", spsa_f_minus))
        cat(sprintf("    gradient: [%.4f, %.4f, %.4f, %.4f, %.4f]\n", 
                    g_hat[1], g_hat[2], g_hat[3], g_hat[4], g_hat[5]))
        cat(sprintf("    Updated theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))
        
        spsa_params$f_plus_history = c(spsa_params$f_plus_history, spsa_f_plus)
        spsa_params$f_minus_history = c(spsa_params$f_minus_history, spsa_f_minus)
        spsa_params$score_history = c(spsa_params$score_history, (spsa_f_plus + spsa_f_minus) / 2)
        spsa_params$theta_history = rbind(spsa_params$theta_history, spsa_params$theta)
        spsa_params$t_history = c(spsa_params$t_history, t)
        
        cat(sprintf("    → Resetting simulation state for next iteration\n"))
        reset_simulation_state()
        spsa_score_window_e = numeric(0)
        spsa_score_window_p = numeric(0)
        
        spsa_delta = sample(c(-1, 1), length(spsa_params$theta), replace = TRUE)
        spsa_params$k = spsa_params$k + 1
        
        c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
        perturb = c_k*spsa_params$c_scale*spsa_delta
        
        theta_plus = clip_theta(
          spsa_params$theta + perturb,
          spsa_params$lower, spsa_params$upper
        )
        
        diffusion_speed_SAMPs          = theta_plus[1]
        add_SAMPs                      = theta_plus[2]
        SAMPs_decay                    = theta_plus[3]
        treg_discrimination_efficiency = theta_plus[4]
        activation_threshold_SAMPs     = theta_plus[5]
        
        cat(sprintf("    → Starting PLUS phase (iter %d) with theta_plus: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    spsa_params$k, theta_plus[1], theta_plus[2], theta_plus[3], theta_plus[4], theta_plus[5]))
        cat(sprintf("--- END SPSA UPDATE ---\n\n"))
        
        spsa_phase = "plus"
      }
    }
  }
  
  # ========================================================================
  # SPSA PARAMETER UPDATE (at regular intervals) - SKIP IF EARLY TRIGGERED
  # ========================================================================
  if (!early_triggered && t > 1 && (t %% spsa_update_interval) == 0) {
    
    current_window_mean_e = mean(spsa_score_window_e)
    current_window_min_e  = min(spsa_score_window_e)
    
    current_window_mean_p = mean(spsa_score_window_p)
    current_window_min_p  = min(spsa_score_window_p)
    
    # objective = 0*current_window_mean + 1*current_window_min ### until try 41
    objective = current_window_min_e + current_window_min_p 
    
    cat(sprintf("\n>>> [REGULAR INTERVAL UPDATE] <<<\n"))
    cat(sprintf("    Time: t=%d | Phase: %s | SPSA iter: %d\n", t, spsa_phase, spsa_params$k))
    cat(sprintf("    Window stats - Epithelial (mean: %.1f, min: %.1f) | Pathogens (mean: %.1f, min: %.1f)\n",
                current_window_mean_e, current_window_min_e, current_window_mean_p, current_window_min_p))
    cat(sprintf("    Objective: %.2f\n", objective))
    
    spsa_score_window_e = numeric(0)
    spsa_score_window_p = numeric(0)
    
    if (spsa_phase == "baseline") {
      spsa_last_window_mean = objective
      spsa_delta = sample(c(-1, 1), length(spsa_params$theta), replace = TRUE)
      spsa_params$k = spsa_params$k + 1
      
      cat(sprintf("    Baseline objective recorded: %.2f\n", objective))
      
      c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
      perturb = c_k*spsa_params$c_scale*spsa_delta
      
      theta_plus = clip_theta(
        spsa_params$theta + perturb,
        spsa_params$lower, spsa_params$upper
      )
      
      diffusion_speed_SAMPs          = theta_plus[1]
      add_SAMPs                      = theta_plus[2]
      SAMPs_decay                    = theta_plus[3]
      treg_discrimination_efficiency = theta_plus[4]
      activation_threshold_SAMPs     = theta_plus[5]
      
      cat(sprintf("    → Switching to PLUS phase with theta_plus: [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n",
                  theta_plus[1], theta_plus[2], theta_plus[3], theta_plus[4], theta_plus[5]))
      
      spsa_phase = "plus"
      
    } else if (spsa_phase == "plus") {
      spsa_f_plus = objective
      cat(sprintf("    Recorded: f_plus = %.2f\n", spsa_f_plus))
      
      cat(sprintf("    → Resetting simulation state for MINUS evaluation\n"))
      reset_simulation_state()
      
      c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
      perturb = c_k*spsa_params$c_scale*spsa_delta
      
      theta_minus = clip_theta(
        spsa_params$theta - perturb,
        spsa_params$lower, spsa_params$upper
      )
      
      diffusion_speed_SAMPs          = theta_minus[1]
      add_SAMPs                      = theta_minus[2]
      SAMPs_decay                    = theta_minus[3]
      treg_discrimination_efficiency = theta_minus[4]
      activation_threshold_SAMPs     = theta_minus[5]
      
      cat(sprintf("    → Switching to MINUS phase with theta_minus: [%.4f, %.4f, %.4f, %.4f, %.4f]\n\n",
                  theta_minus[1], theta_minus[2], theta_minus[3], theta_minus[4], theta_minus[5]))
      
      spsa_phase = "minus"
      
    } else if (spsa_phase == "minus") {
      spsa_f_minus = objective
      cat(sprintf("    Recorded: f_minus = %.2f\n", spsa_f_minus))
      
      c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
      a_k = spsa_params$a / (spsa_params$A + spsa_params$k)^spsa_params$alpha
      
      perturb = c_k*spsa_params$c_scale*spsa_delta
      g_hat = (spsa_f_plus - spsa_f_minus) / (2*perturb)
      
      step = a_k*spsa_params$a_scale*g_hat
      
      spsa_params$theta = clip_theta(
        spsa_params$theta + step,
        spsa_params$lower, spsa_params$upper
      )
      
      diffusion_speed_SAMPs          = spsa_params$theta[1]
      add_SAMPs                      = spsa_params$theta[2]
      SAMPs_decay                    = spsa_params$theta[3]
      treg_discrimination_efficiency = spsa_params$theta[4]
      activation_threshold_SAMPs     = spsa_params$theta[5]
      
      cat(sprintf("\n--- SPSA UPDATE (REGULAR, iter %d) ---\n", spsa_params$k))
      cat(sprintf("    f_plus:  %.2f\n", spsa_f_plus))
      cat(sprintf("    f_minus: %.2f\n", spsa_f_minus))
      cat(sprintf("    gradient: [%.4f, %.4f, %.4f, %.4f, %.4f]\n", 
                  g_hat[1], g_hat[2], g_hat[3], g_hat[4], g_hat[5]))
      cat(sprintf("    Updated theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                  spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                  spsa_params$theta[4], spsa_params$theta[5]))
      
      spsa_params$f_plus_history = c(spsa_params$f_plus_history, spsa_f_plus)
      spsa_params$f_minus_history = c(spsa_params$f_minus_history, spsa_f_minus)
      spsa_params$score_history = c(spsa_params$score_history, (spsa_f_plus + spsa_f_minus) / 2)
      spsa_params$theta_history = rbind(spsa_params$theta_history, spsa_params$theta)
      spsa_params$t_history = c(spsa_params$t_history, t)
      
      cat(sprintf("    → Resetting simulation state for next iteration\n"))
      reset_simulation_state()
      
      spsa_delta = sample(c(-1, 1), length(spsa_params$theta), replace = TRUE)
      spsa_params$k = spsa_params$k + 1
      
      c_k = spsa_params$c / (spsa_params$k)^spsa_params$gamma
      perturb = c_k*spsa_params$c_scale*spsa_delta
      
      theta_plus = clip_theta(
        spsa_params$theta + perturb,
        spsa_params$lower, spsa_params$upper
      )
      
      diffusion_speed_SAMPs          = theta_plus[1]
      add_SAMPs                      = theta_plus[2]
      SAMPs_decay                    = theta_plus[3]
      treg_discrimination_efficiency = theta_plus[4]
      activation_threshold_SAMPs     = theta_plus[5]
      
      cat(sprintf("    → Starting PLUS phase (iter %d) with theta_plus: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                  spsa_params$k, theta_plus[1], theta_plus[2], theta_plus[3], theta_plus[4], theta_plus[5]))
      cat(sprintf("--- END SPSA UPDATE ---\n\n"))
      
      spsa_phase = "plus"
    }
  }
  source("./MISC/MAIN_SIM.R")
}

# ============================================================================
# SAVE LONGITUDINAL DATA
# ============================================================================
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

# ============================================================================
# AFTER SIMULATION: Compile SPSA optimization history
# ============================================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("SPSA OPTIMIZATION COMPLETE\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Total iterations: %d\n", spsa_params$k))
cat(sprintf("Final theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
            spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
            spsa_params$theta[4], spsa_params$theta[5]))
cat(sprintf("========================================\n\n"))

spsa_results = data.frame(
  iter = seq_along(spsa_params$score_history),
  t = spsa_params$t_history,
  f_plus = spsa_params$f_plus_history,
  f_minus = spsa_params$f_minus_history,
  score_mean = spsa_params$score_history,
  diffusion_speed_SAMPs = spsa_params$theta_history[, 1],
  add_SAMPs = spsa_params$theta_history[, 2],
  SAMPs_decay = spsa_params$theta_history[, 3],
  treg_discrimination_efficiency = spsa_params$theta_history[, 4],
  activation_threshold_SAMPs = spsa_params$theta_history[, 5]
)

# write.csv(spsa_results, 
#           file = paste0("./spsa_optimization_", param_set_id_use, ".csv"),
#           row.names = FALSE)
