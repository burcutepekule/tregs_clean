# ============================================================================
# ODE VERSION - Replaces RUN_REPS_CPP_ABM_PAMPS.R
# ============================================================================
# This script runs ODE simulations analogous to the ABM version

source('./MISC/ODE_SYSTEM.R')

# ========================================================================
# PREPARE ODE PARAMETERS
# ========================================================================
# Create parameter list from assigned variables
# Note: Some ABM-specific parameters need to be adapted for ODE use

ode_params = list(
  # === Grid/environment parameters ===
  grid_size         = grid_size,
  E_max             = 1,  # Maximum epithelial health score - let's assume max is 1 for the sake of making it east to interpret?
  n_phagocytes      = n_phagocytes,
  n_commensals_lp   = n_commensals_lp,
  injury_percentage = injury_percentage,

  # === Macrophage polarization ===
  activation_threshold_danger = activation_threshold_danger,
  activation_threshold_SAMPs  = activation_threshold_SAMPs,
  polarization_speed = .1,  # ODE-specific: rate of M0->M1/M2 conversion

  # === Macrophage recruitment ===
  recruitment_rate_danger = recruitment_rate_danger,
  
  # === Macrophage activities ===
  activity_ROS_M0_baseline = activity_ROS_M0_baseline,
  activity_ROS_M1_baseline = activity_ROS_M1_baseline,
  activity_ROS_M2_baseline = activity_ROS_M2_baseline,
  activity_engulf_M0_baseline = activity_engulf_M0_baseline,
  activity_engulf_M1_baseline = activity_engulf_M1_baseline,
  activity_engulf_M2_baseline = activity_engulf_M2_baseline,
  
  # === Epithelial dynamics ===
  epith_recovery_chance  = epith_recovery_chance,
  damage_rate_pathogen   = 0.5,  # ODE-specific: pathogen damage rate
  damage_rate_ROS        = 0.1,  # ODE-specific: ROS damage rate
  th_ROS_epith_injury    = th_ROS_epith_injury,
  
  # === Microbe leakage ===
  rate_leak_pathogen_injury    = rate_leak_pathogen_injury,
  rate_leak_commensal_baseline = rate_leak_commensal_baseline,
  rate_leak_commensal_injury   = rate_leak_commensal_injury,
  
  # === Microbe killing ===
  th_ROS_microbe = th_ROS_microbe,
  kill_rate_ROS  = 1.0,  # ODE-specific: ROS killing efficiency
  K_phagocytosis = 100,  # ODE-specific: half-saturation for phagocytosis
  
  # === Treg parameters ===
  allow_tregs                 = allow_tregs,
  treg_activation_efficiency  = 0.01,  # ODE-specific: Treg activation rate
  K_treg_activation           = 50,  # ODE-specific: half-saturation for Treg activation
  
  # === Signaling molecule production ===
  add_ROS   = add_ROS,
  add_DAMPs = add_DAMPs,
  add_PAMPs = add_PAMPs,
  add_SAMPs = add_SAMPs,
  
  # === Signaling molecule decay ===
  ros_decay   = ros_decay,
  DAMPs_decay = DAMPs_decay,
  PAMPs_decay = PAMPs_decay,
  SAMPs_decay = SAMPs_decay,

  # === Epithelial recovery constants ==
  k_recovery = 100
  
)

# ========================================================================
# RUN ODE SIMULATION
# ========================================================================
cat(sprintf("Running ODE simulation: param_set=%d, pat_level=%d\n",
            param_set_use$param_set_id, pat_level))

# Solve ODEs
ode_output = run_ode_simulation(
  parameters = ode_params,
  t_max = t_max,
  dt = 1  # Time step (same as ABM)
)

# ========================================================================
# PROCESS OUTPUT
# ========================================================================
scenario_info = list(
  param_set_id = param_set_use$param_set_id,
  sterile = sterile,
  macspec_on = macspec_on,
  allow_tregs = allow_tregs,
  ros_level = ros_level,
  pat_level = pat_level,
  randomize_tregs = randomize_tregs
)

longitudinal_df = process_ode_output(ode_output, ode_params, scenario_info)

# ========================================================================
# APPEND TO MASTER DATAFRAME
# ========================================================================
longitudinal_df_keep = rbind(longitudinal_df_keep, longitudinal_df)

cat(sprintf("  Completed. Final epithelial score: %.2f, Final pathogens: %.2f\n",
            tail(longitudinal_df$epithelial_score, 1),
            tail(longitudinal_df$pathogen, 1)))

