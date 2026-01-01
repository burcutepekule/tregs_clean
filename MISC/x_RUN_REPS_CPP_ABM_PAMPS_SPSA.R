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
  upper = c(0.120, 0.500, 0.500, 1.000, 1.000),
  
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

# Detection settings
collapse_threshold = 30
collapse_duration  = 25
collapse_rate      = 0.75

### from try 25 to 40
# success_threshold  = 120 # used to be 100
# success_duration   = 125 # used to be 30
# success_rate       = 0.95

### try 41
# success_threshold  = 145 # used to be 100
# success_duration   = 250 # used to be 30
# success_rate       = 0.95

success_threshold_e= 150*0.75 
success_threshold_p= 10
success_duration   = 250 
success_rate       = 0.95

# Initialize success log file
cat("# SPSA Success Log\n", file = paste0("./spsa_successes_", param_set_id_use,"_scenario_",scenario_ind,".txt"), append = FALSE)

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
      pct_below_threshold = sum(spsa_score_window_e < collapse_threshold) / length(spsa_score_window_e)
      is_collapse = all(recent_scores_collapse < collapse_threshold) || 
        (length(spsa_score_window_e) >= collapse_duration && pct_below_threshold > collapse_rate)
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
        early_objective = mean(spsa_score_window_e)
        cat(sprintf("  [COLLAPSE at t=%d, score=%.1f, phase=%s]\n", 
                    t, current_epithelial_score, spsa_phase))
      } else {
        early_objective = mean(spsa_score_window_e)
        cat(sprintf("  [SUCCESS at t=%d, score=%.1f, phase=%s]\n", 
                    t, current_epithelial_score, spsa_phase))
        cat(sprintf("  theta: [%.3f, %.3f, %.3f, %.3f, %.3f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))
        
        # Write success to file
        success_line = sprintf("t=%d | iter=%d | score=%.1f | phase=%s | theta=[%.4f, %.4f, %.4f, %.4f, %.4f]\n",
                               t, spsa_params$k, early_objective, spsa_phase,
                               spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                               spsa_params$theta[4], spsa_params$theta[5])
        
        cat(success_line, file = paste0("./spsa_successes_", param_set_id_use,"_scenario_",scenario_ind,".txt"), append = TRUE)
        
        # Continue optimization like a collapse (don't break)
      }
      
      if (spsa_phase == "plus") {
        spsa_f_plus = early_objective
        
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
        
        spsa_phase = "minus"
        
      } else if (spsa_phase == "minus") {
        spsa_f_minus = early_objective
        
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
        
        cat(sprintf("t=%d | SPSA iter %d | f+: %.1f, f-: %.1f [EARLY]\n",
                    t, spsa_params$k, spsa_f_plus, spsa_f_minus))
        cat(sprintf("  theta: [%.3f, %.3f, %.3f, %.3f, %.3f]\n",
                    spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                    spsa_params$theta[4], spsa_params$theta[5]))
        
        spsa_params$f_plus_history = c(spsa_params$f_plus_history, spsa_f_plus)
        spsa_params$f_minus_history = c(spsa_params$f_minus_history, spsa_f_minus)
        spsa_params$score_history = c(spsa_params$score_history, (spsa_f_plus + spsa_f_minus) / 2)
        spsa_params$theta_history = rbind(spsa_params$theta_history, spsa_params$theta)
        spsa_params$t_history = c(spsa_params$t_history, t)
        
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
    
    spsa_score_window_e = numeric(0)
    spsa_score_window_p = numeric(0)
    
    if (spsa_phase == "baseline") {
      spsa_last_window_mean = objective
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
      
      spsa_phase = "plus"
      
    } else if (spsa_phase == "plus") {
      spsa_f_plus = objective
      
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
      
      spsa_phase = "minus"
      
    } else if (spsa_phase == "minus") {
      spsa_f_minus = objective
      
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
      
      cat(sprintf("t=%d | SPSA iter %d | f+: %.1f, f-: %.1f\n",
                  t, spsa_params$k, spsa_f_plus, spsa_f_minus))
      cat(sprintf("  theta: [%.3f, %.3f, %.3f, %.3f, %.3f]\n",
                  spsa_params$theta[1], spsa_params$theta[2], spsa_params$theta[3],
                  spsa_params$theta[4], spsa_params$theta[5]))
      
      spsa_params$f_plus_history = c(spsa_params$f_plus_history, spsa_f_plus)
      spsa_params$f_minus_history = c(spsa_params$f_minus_history, spsa_f_minus)
      spsa_params$score_history = c(spsa_params$score_history, (spsa_f_plus + spsa_f_minus) / 2)
      spsa_params$theta_history = rbind(spsa_params$theta_history, spsa_params$theta)
      spsa_params$t_history = c(spsa_params$t_history, t)
      
      reset_simulation_state()
      cat("  [State reset]\n")
      
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
      
      spsa_phase = "plus"
    }
  }
  
  # Update injury site
  injury_site_updated = which(epithelium$level_injury > 0)
  
  # ========================================================================
  # UPDATE SAMPs (from activated Tregs)
  # ========================================================================
  active_tregs = which(treg_phenotype == 1)
  if (length(active_tregs) > 0) {
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
  DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury*add_DAMPs
  DAMPs[1, ] = DAMPs[1, ] + logistic_scaled_0_to_5_quantized(
    pathogen_epithelium_counts + commensal_epithelium_counts
  )*add_DAMPs
  
  # ========================================================================
  # UPDATE PAMPs 
  # ========================================================================
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
  # ========================================================================
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
    
    for (x in 1:grid_size) {
      danger_at_border = danger_signal_grid[grid_size, x]
      n_recruit = rpois(1, lambda = recruitment_rate_danger*danger_at_border)
      
      if (n_recruit > 0) {
        phagocyte_x = c(phagocyte_x, rep(x, n_recruit))
        phagocyte_y = c(phagocyte_y, rep(grid_size, n_recruit))
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
  density_matrix_phagocytes = DAMPs + PAMPs
  
  all_equal_treg = all(density_matrix_tregs == density_matrix_tregs[1, 1])
  all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])
  
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
  # UPDATE PHAGOCYTE PHENOTYPES
  # ========================================================================
  M0_indices = which(phagocyte_phenotype == 0)
  M1_indices = which(phagocyte_phenotype == 1)
  M2_indices = which(phagocyte_phenotype == 2)
  
  phagocyte_commensals_engulfed = rowSums(phagocyte_bacteria_registry > 0)
  phagocyte_pathogens_engulfed = rowSums(phagocyte_bacteria_registry < 0)
  
  if (length(M0_indices) > 0) {
    signals = calculate_phagocyte_signals_cpp(
      M0_indices, phagocyte_x, phagocyte_y,
      act_radius_DAMPs, act_radius_SAMPs, DAMPs, SAMPs, grid_size
    )
    avg_DAMPs_vec = signals$avg_DAMPs
    avg_SAMPs_vec = signals$avg_SAMPs
    
    pamps_signals = calculate_phagocyte_signals_cpp(
      M0_indices, phagocyte_x, phagocyte_y,
      act_radius_DAMPs, act_radius_DAMPs, PAMPs, PAMPs, grid_size
    )
    avg_PAMPs_vec = pamps_signals$avg_DAMPs
    
    for (idx in seq_along(M0_indices)) {
      i = M0_indices[idx]
      avg_DAMPs = avg_DAMPs_vec[idx]
      avg_SAMPs = avg_SAMPs_vec[idx]
      avg_PAMPs = avg_PAMPs_vec[idx]
      
      danger_signal = avg_DAMPs + avg_PAMPs
      
      if (danger_signal >= activation_threshold_danger && danger_signal > avg_SAMPs) {
        phagocyte_phenotype[i] = 1
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M1_baseline
        phagocyte_activity_engulf[i] = activity_engulf_M1_baseline
      } 
      ## Prevent M0 → M2 transitions so that Tregs can only suppress already-activated M1 macrophages (M1 → M2), not activate naive M0 macrophages.
      # else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > danger_signal) {
      #   phagocyte_phenotype[i] = 2
      #   phagocyte_active_age[i] = 1
      #   phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
      #   phagocyte_activity_engulf[i] = activity_engulf_M2_baseline
      # }
    }
  }
  
  active_indices = c(M1_indices, M2_indices)
  if (length(active_indices) > 0) {
    phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
    old_enough = phagocyte_active_age[active_indices] >= active_age_limit
    candidates = active_indices[old_enough]
    
    if (length(candidates) > 0) {
      signals = calculate_phagocyte_signals_cpp(
        candidates, phagocyte_x, phagocyte_y,
        act_radius_DAMPs, act_radius_SAMPs, DAMPs, SAMPs, grid_size
      )
      avg_DAMPs_vec = signals$avg_DAMPs
      avg_SAMPs_vec = signals$avg_SAMPs
      
      pamps_signals = calculate_phagocyte_signals_cpp(
        candidates, phagocyte_x, phagocyte_y,
        act_radius_DAMPs, act_radius_DAMPs, PAMPs, PAMPs, grid_size
      )
      avg_PAMPs_vec = pamps_signals$avg_DAMPs
      
      for (idx in seq_along(candidates)) {
        i = candidates[idx]
        avg_DAMPs = avg_DAMPs_vec[idx]
        avg_SAMPs = avg_SAMPs_vec[idx]
        avg_PAMPs = avg_PAMPs_vec[idx]
        
        danger_signal = avg_DAMPs + avg_PAMPs
        
        if(macspec_on > 0){
          num_pat_engulfed = phagocyte_pathogens_engulfed[i]
          num_com_engulfed = phagocyte_commensals_engulfed[i]
          
          DAMPs_dominant = (danger_signal >= activation_threshold_danger && danger_signal > avg_SAMPs)
          SAMPs_dominant = (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > danger_signal)
          
          pathogen_engulfment_dominant  = FALSE
          commensal_engulfment_dominant = FALSE
          
          if ((num_pat_engulfed + num_com_engulfed) > 0) {
            rat_com_pat_real = num_com_engulfed / (num_com_engulfed + num_pat_engulfed)
            
            commensal_presented = runif(1) < rat_com_pat_real
            
            if (commensal_presented) {
              mac_identifies_as_commensal = runif(1) < mac_discrimination_efficiency
            } else {
              mac_identifies_as_commensal = runif(1) < (1 - mac_discrimination_efficiency)
            }
            
            pathogen_engulfment_dominant  = !mac_identifies_as_commensal
            commensal_engulfment_dominant = mac_identifies_as_commensal
          }

          # FIX: Phenotype-specific deactivation logic (macspec mode)
          current_phenotype = phagocyte_phenotype[i]

          if (avg_SAMPs < activation_threshold_SAMPs && danger_signal < activation_threshold_danger) {
            # Both signals low → always deactivate
            phagocyte_phenotype[i] = 0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          } else if (current_phenotype == 1 && danger_signal < activation_threshold_danger) {
            # M1 deactivates when danger drops, regardless of SAMPs
            phagocyte_phenotype[i] = 0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          } else if (DAMPs_dominant || pathogen_engulfment_dominant) {
            phagocyte_phenotype[i] = 1
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M1_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M1_baseline
          } else if (SAMPs_dominant && commensal_engulfment_dominant) {
            phagocyte_phenotype[i] = 2
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M2_baseline
          }
        } else {
          # FIX: Phenotype-specific deactivation logic (vanilla mode)
          current_phenotype = phagocyte_phenotype[i]

          if (avg_SAMPs < activation_threshold_SAMPs && danger_signal < activation_threshold_danger) {
            # Both signals low → always deactivate
            phagocyte_phenotype[i] = 0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          } else if (current_phenotype == 1 && danger_signal < activation_threshold_danger) {
            # M1 deactivates when danger drops, regardless of SAMPs
            phagocyte_phenotype[i] = 0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          } else if (danger_signal >= activation_threshold_danger && danger_signal > avg_SAMPs) {
            # Danger dominant → M1
            phagocyte_phenotype[i] = 1
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M1_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M1_baseline
          } else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > danger_signal) {
            # SAMPs dominant → M2
            phagocyte_phenotype[i] = 2
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M2_baseline
          }
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
    
    if (nrow(pathogen_coords) > 0) {
      pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
      pathogen_indices = which(pathogen_overlap)
      
      if (length(pathogen_indices) > 0) {
        engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
        indices_to_engulf = pathogen_indices[engulf_success]
        
        if (length(indices_to_engulf) > 0) {
          pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
          
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
    
    if (nrow(commensal_coords) > 0) {
      commensal_overlap = (commensal_coords[, "x"] == px) & (commensal_coords[, "y"] == py)
      commensal_indices = which(commensal_overlap)
      
      if (length(commensal_indices) > 0) {
        engulf_success = runif(length(commensal_indices)) < phagocyte_activity_engulf[i]
        indices_to_engulf = commensal_indices[engulf_success]
        
        if (length(indices_to_engulf) > 0) {
          commensal_coords = commensal_coords[-indices_to_engulf, , drop = FALSE]
          
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
      phagocyte_bacteria_registry[i, ] = shift_insert_fast_cpp(
        phagocyte_bacteria_registry[i, ],
        numeric(0)
      )
    }
  }
  
  # ========================================================================
  # CALCULATE COUNTS FROM REGISTRY
  # ========================================================================
  phagocyte_commensals_engulfed = rowSums(phagocyte_bacteria_registry > 0)
  phagocyte_pathogens_engulfed  = rowSums(phagocyte_bacteria_registry < 0)
  
  # ========================================================================
  # TREG ACTIVATION & EFFECTOR ACTIONS
  # ========================================================================
  M1_phagocyte_indices = which(phagocyte_phenotype == 1)
  M2_phagocyte_indices = which(phagocyte_phenotype == 2)
  M_activate_phagocyte_indices = c(M1_phagocyte_indices, M2_phagocyte_indices)
  
  if (allow_tregs == 1 && length(M_activate_phagocyte_indices) > 0) {
    for (i in M_activate_phagocyte_indices) {
      px = phagocyte_x[i]
      py = phagocyte_y[i]
      
      nearby_treg_indices = find_nearby_tregs_cpp(
        px, py, treg_x, treg_y, treg_vicinity_effect
      )
      
      if (length(nearby_treg_indices) > 0) {
        num_pat_antigens = phagocyte_pathogens_engulfed[i]
        num_com_antigens = phagocyte_commensals_engulfed[i]
        
        if ((num_pat_antigens + num_com_antigens) > 0) {
          rat_com_pat_real = num_com_antigens / (num_com_antigens + num_pat_antigens)
          
          commensal_presented = runif(1) < rat_com_pat_real
          
          if (commensal_presented) {
            treg_identifies_as_commensal = runif(1) < treg_discrimination_efficiency
          } else {
            treg_identifies_as_commensal = runif(1) < (1 - treg_discrimination_efficiency)
          }
          
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
  
  epithelium$level_injury = epithelium$level_injury +
    logistic_scaled_0_to_5_quantized(pathogen_epithelium_counts)
  
  epithelium$level_injury = epithelium$level_injury + as.integer(ros_means > th_ROS_epith_injury)
  
  epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)
  
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

write.csv(spsa_results, 
          file = paste0("./spsa_optimization_", param_set_id_use, ".csv"),
          row.names = FALSE)