# ODE System for Immune Response Model
# Translates the spatial ABM into a well-mixed ODE system
#
# State variables:
#   M0, M1, M2: Macrophage phenotypes
#   P: Pathogen population
#   C: Commensal population
#   E: Epithelial health score (0 = completely damaged, E_max = fully healthy)
#   ROS: Reactive oxygen species concentration
#   DAMPS: Damage-associated molecular patterns
#   PAMPS: Pathogen-associated molecular patterns
#   SAMPS: Suppression-associated molecular patterns

library(deSolve)

# ODE system function
# state: named vector of state variables
# time: current time
# parameters: list of model parameters
ode_system <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    # === Calculate signals for macrophage polarization ===
    danger_signal <- DAMPS + PAMPS
    safety_signal <- SAMPS

    # === Macrophage polarization rates ===
    # M0 polarizes to M1 when danger > safety and exceeds threshold
    polarization_rate_M1 <- 0
    if (danger_signal >= activation_threshold_danger && danger_signal > safety_signal && M0 > 0) {
      # Polarization rate proportional to M0 and signal strength
      polarization_rate_M1 <- polarization_speed * M0 * (danger_signal / (danger_signal + K_polarization))
    }

    # M0 polarizes to M2 when safety > danger and exceeds threshold
    polarization_rate_M2 <- 0
    if (safety_signal >= activation_threshold_SAMPs && safety_signal > danger_signal && M0 > 0) {
      polarization_rate_M2 <- polarization_speed * M0 * (safety_signal / (safety_signal + K_polarization))
    }

    # Macrophage recruitment (danger-dependent)
    recruitment_rate <- recruitment_rate_danger * danger_signal

    # === Epithelial health dynamics ===
    # E ranges from 0 (fully damaged) to E_max (fully healthy)
    # Injury is inverse of health
    injury_fraction <- 1 - (E / E_max)  # 0 = healthy, 1 = fully damaged

    # Epithelial recovery (proportional to current damage)
    recovery_rate <- epith_recovery_chance * (E_max - E)

    # Epithelial damage from pathogens (saturating function)
    damage_from_pathogens <- damage_rate_pathogen * P * E / (K_pathogen_damage + P)

    # Epithelial damage from ROS (above threshold)
    damage_from_ROS <- 0
    if (ROS > th_ROS_epith_injury) {
      damage_from_ROS <- damage_rate_ROS * (ROS - th_ROS_epith_injury) * E / E_max
    }

    # === Microbe dynamics ===
    # Pathogen leakage (injury-dependent)
    pathogen_leakage <- rate_leak_pathogen_injury * injury_fraction * grid_size

    # Commensal leakage (baseline + injury-enhanced)
    commensal_leakage <- rate_leak_commensal_baseline * grid_size +
                         rate_leak_commensal_injury * injury_fraction * grid_size

    # Microbe killing by ROS (above threshold)
    pathogen_kill_ROS <- 0
    commensal_kill_ROS <- 0
    if (ROS > th_ROS_microbe) {
      # Kill rate proportional to ROS above threshold
      ros_kill_efficiency <- kill_rate_ROS * (ROS - th_ROS_microbe)
      pathogen_kill_ROS <- ros_kill_efficiency * P
      commensal_kill_ROS <- ros_kill_efficiency * C
    }

    # Microbe killing by phagocytosis
    # Total phagocytosis capacity
    phagocytosis_M0 <- activity_engulf_M0_baseline * M0
    phagocytosis_M1 <- activity_engulf_M1_baseline * M1
    phagocytosis_M2 <- activity_engulf_M2_baseline * M2
    total_phagocytosis <- phagocytosis_M0 + phagocytosis_M1 + phagocytosis_M2

    # Distribute phagocytosis between pathogens and commensals
    # Based on relative abundances and encounter rates
    total_microbes <- P + C
    if (total_microbes > 0) {
      pathogen_fraction <- P / total_microbes
      commensal_fraction <- C / total_microbes
    } else {
      pathogen_fraction <- 0
      commensal_fraction <- 0
    }

    pathogen_kill_phago <- total_phagocytosis * pathogen_fraction * P / (K_phagocytosis + total_microbes)
    commensal_kill_phago <- total_phagocytosis * commensal_fraction * C / (K_phagocytosis + total_microbes)

    # === Treg activation and SAMP production ===
    # Tregs activate when macrophages have engulfed bacteria
    # Simplified: Treg activation proportional to M1+M2 and commensal abundance
    if (allow_tregs == 1) {
      treg_activation_rate <- treg_activation_efficiency * (M1 + M2) * C / (K_treg_activation + C)
    } else {
      treg_activation_rate <- 0
    }

    # === Signaling molecule dynamics ===
    # ROS production by M1 macrophages
    ROS_production <- activity_ROS_M1_baseline * M1 * add_ROS
    ROS_decay_rate <- ros_decay * ROS

    # DAMPS production from injured epithelium and microbes
    DAMPS_production_injury <- injury_fraction * add_DAMPs * grid_size
    # Microbe-induced DAMPs (saturating)
    DAMPS_production_microbes <- add_DAMPs * (P + C) / (K_DAMPS_microbe + P + C)
    DAMPS_decay_rate <- DAMPs_decay * DAMPS

    # PAMPS production by pathogens
    PAMPS_production <- add_PAMPs * P
    PAMPS_decay_rate <- PAMPs_decay * PAMPS

    # SAMPS production by activated Tregs
    SAMPS_production <- add_SAMPs * treg_activation_rate
    SAMPS_decay_rate <- SAMPs_decay * SAMPS

    # === Differential equations ===
    dM0 <- recruitment_rate - polarization_rate_M1 - polarization_rate_M2
    dM1 <- polarization_rate_M1
    dM2 <- polarization_rate_M2

    dP <- pathogen_leakage - pathogen_kill_ROS - pathogen_kill_phago
    dC <- commensal_leakage - commensal_kill_ROS - commensal_kill_phago

    dE <- recovery_rate - damage_from_pathogens - damage_from_ROS

    dROS <- ROS_production - ROS_decay_rate
    dDAMPS <- DAMPS_production_injury + DAMPS_production_microbes - DAMPS_decay_rate
    dPAMPS <- PAMPS_production - PAMPS_decay_rate
    dSAMPS <- SAMPS_production - SAMPS_decay_rate

    # Return derivatives
    list(c(dM0, dM1, dM2, dP, dC, dE, dROS, dDAMPS, dPAMPS, dSAMPS),
         # Also return some useful metrics for output
         danger_signal = danger_signal,
         safety_signal = safety_signal,
         injury_fraction = injury_fraction,
         pathogen_kill_ROS = pathogen_kill_ROS,
         commensal_kill_ROS = commensal_kill_ROS,
         pathogen_kill_phago = pathogen_kill_phago,
         commensal_kill_phago = commensal_kill_phago)
  })
}

# Function to set up initial conditions
setup_initial_conditions <- function(parameters) {
  with(as.list(parameters), {
    # Initial macrophage population (all M0)
    M0_init <- n_phagocytes
    M1_init <- 0
    M2_init <- 0

    # Initial microbes in lamina propria
    P_init <- 0  # Start with no pathogens (they leak in)
    C_init <- n_commensals_lp  # Initial commensals

    # Initial epithelial health
    # Start with injury_percentage of epithelium damaged
    # E_max represents fully healthy epithelium
    E_max <- grid_size * grid_size * 6  # Max score (all cells at health level 6)
    injury_level_initial <- max_level_injury * (injury_percentage / 100)
    E_init <- E_max * (1 - injury_level_initial / max_level_injury)

    # Initial signaling molecules (start at zero)
    ROS_init <- 0
    DAMPS_init <- 0
    PAMPS_init <- 0
    SAMPS_init <- 0

    state <- c(M0 = M0_init, M1 = M1_init, M2 = M2_init,
               P = P_init, C = C_init, E = E_init,
               ROS = ROS_init, DAMPS = DAMPS_init,
               PAMPS = PAMPS_init, SAMPS = SAMPS_init)

    return(state)
  })
}

# Function to run ODE simulation
run_ode_simulation <- function(parameters, t_max = 5000, dt = 1) {
  # Set up initial conditions
  initial_state <- setup_initial_conditions(parameters)

  # Time points
  times <- seq(0, t_max, by = dt)

  # Solve ODE system
  out <- ode(y = initial_state,
             times = times,
             func = ode_system,
             parms = parameters,
             method = "lsoda")  # Adaptive solver for stiff systems

  # Convert to data frame
  out_df <- as.data.frame(out)

  return(out_df)
}

# Function to process ODE output into ABM-like format
# This creates outputs similar to the ABM for consistency
process_ode_output <- function(ode_output, parameters, scenario_info) {
  df <- ode_output

  # Rename time column
  names(df)[names(df) == "time"] <- "t"

  # Calculate epithelial score (for compatibility with ABM output)
  # E is already the epithelial score
  df$epithelial_score <- df$E

  # Calculate injury fraction for reference
  df$injury_fraction <- 1 - (df$E / (parameters$grid_size^2 * 6))

  # Add cumulative death counters (integrate killing rates over time)
  df$C_ROS <- cumsum(df$commensal_kill_ROS)
  df$C_M0 <- 0  # Not tracking separately in ODE
  df$C_M1 <- cumsum(df$commensal_kill_phago * df$M1 / (df$M0 + df$M1 + df$M2 + 1e-10))
  df$C_M2 <- cumsum(df$commensal_kill_phago * df$M2 / (df$M0 + df$M1 + df$M2 + 1e-10))

  df$P_ROS <- cumsum(df$pathogen_kill_ROS)
  df$P_M0 <- 0
  df$P_M1 <- cumsum(df$pathogen_kill_phago * df$M1 / (df$M0 + df$M1 + df$M2 + 1e-10))
  df$P_M2 <- cumsum(df$pathogen_kill_phago * df$M2 / (df$M0 + df$M1 + df$M2 + 1e-10))

  # Rename compartments to match ABM output
  df$phagocyte_M0 <- df$M0
  df$phagocyte_M1 <- df$M1
  df$phagocyte_M2 <- df$M2
  df$pathogen <- df$P
  df$commensal <- df$C

  # Add scenario metadata
  df$param_set_id <- scenario_info$param_set_id
  df$rep_id <- scenario_info$rep_id
  df$sterile <- scenario_info$sterile
  df$macspec_on <- scenario_info$macspec_on
  df$tregs_on <- scenario_info$allow_tregs
  df$ros_level <- scenario_info$ros_level
  df$pat_level <- scenario_info$pat_level
  df$randomize_tregs <- scenario_info$randomize_tregs

  # Calculate steady state times (simple approach: when derivative < threshold)
  # For epithelial score
  df$epithelial_score_change <- c(0, diff(df$epithelial_score))
  steady_state_threshold <- 0.01
  ss_e_idx <- which(abs(df$epithelial_score_change) < steady_state_threshold)
  if (length(ss_e_idx) > 0) {
    df$time_ss_e <- df$t[ss_e_idx[1]]
  } else {
    df$time_ss_e <- NA
  }

  # For pathogen
  df$pathogen_change <- c(0, diff(df$pathogen))
  ss_p_idx <- which(abs(df$pathogen_change) < steady_state_threshold)
  if (length(ss_p_idx) > 0) {
    df$time_ss_p <- df$t[ss_p_idx[1]]
  } else {
    df$time_ss_p <- NA
  }

  return(df)
}
