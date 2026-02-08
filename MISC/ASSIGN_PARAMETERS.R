# ============================================================================
# ASSIGN PARAMETERS FROM CSV (WITH MACROPHAGE SPECIFICITY)
# ============================================================================
# Thresholds
th_ROS_microbe        = param_set_use$th_ROS_microbe
th_ROS_epith_injury   = param_set_use$th_ROS_epith_injury
epith_recovery_chance = param_set_use$epith_recovery_chance

# Diffusion speeds
diffusion_speed_DAMPs = param_set_use$diffusion_speed_DAMPs
diffusion_speed_SAMPs = param_set_use$diffusion_speed_SAMPs
diffusion_speed_PAMPs = param_set_use$diffusion_speed_PAMPs
diffusion_speed_ROS   = param_set_use$diffusion_speed_ROS

# ==============================================================================
# ====================== OVERWRITE - easiest for now
add_ROS                   = ros_level*0.1 # ros_level=0, control, until 10, which is max 1=10*0.1
rate_leak_pathogen_injury = ifelse(sterile == 1, 0.0, pat_level)
# ==============================================================================

add_DAMPs = param_set_use$add_DAMPs
add_SAMPs = param_set_use$add_SAMPs
add_PAMPs = param_set_use$add_PAMPs


##### PLAYING WITH THINGS HERE, BE CAREFUL!
# Decay rates
ros_decay   = 10*param_set_use$ros_decay
DAMPs_decay = param_set_use$DAMPs_decay
SAMPs_decay = param_set_use$SAMPs_decay
PAMPs_decay = param_set_use$PAMPs_decay

# Activation thresholds
activation_threshold_danger = param_set_use$activation_threshold_danger
activation_threshold_SAMPs  = param_set_use$activation_threshold_SAMPs
activation_threshold_SAMPs  = 2*activation_threshold_SAMPs

# Engulfment activities
activity_engulf_M0_baseline = param_set_use$activity_engulf_M0_baseline
activity_engulf_M1_baseline = param_set_use$activity_engulf_M1_baseline
activity_engulf_M2_baseline = m2_on*param_set_use$activity_engulf_M2_baseline

# ROS production activities
activity_ROS_M0_baseline = 0
activity_ROS_M2_baseline = 0
activity_ROS_M1_baseline = param_set_use$activity_ROS_M1_baseline

# Leak rates
rate_leak_commensal_injury   = param_set_use$rate_leak_commensal_injury
rate_leak_commensal_baseline = param_set_use$rate_leak_commensal_baseline

# Phagocyte parameters
active_age_limit = as.integer(param_set_use$active_age_limit)
cc_phagocyte     = as.integer(param_set_use$cc_phagocyte)
# digestion_time   = 1 # by default digestion happens every time point

# Treg parameters
treg_vicinity_effect = 1
treg_discrimination_efficiency = param_set_use$treg_discrimination_efficiency

# Macrophage specificity parameters
if(macspec_on == 2){# PERFECT DISCRIMINATION
  mac_discrimination_efficiency = 1 
}else if(macspec_on == 1){# SAME DISCRIMINATION - Fair comparison
  mac_discrimination_efficiency = treg_discrimination_efficiency 
}else{
  mac_discrimination_efficiency = 0 # no need to have this, shouldn't throw an error
}

recruitment_rate_danger   = param_set_use$recruitment_rate_danger

# ============================================================================
# INITIALIZE SIMULATION
# ============================================================================
injury_site = get_middle_percent(seq(1, grid_size), injury_percentage)
n_pathogens_lp = round(rate_leak_pathogen_injury*length(injury_site))

# ============================================================================
# DIFFUSION MODEL PARAMETERS (hardcoded - bacteria/cells slower than signals)
# ============================================================================
# # Microbe diffusion (slower than signaling molecules)
# diffusion_speed_microbe = 0.1   # ~10x slower than typical signal diffusion
# decay_rate_microbe = 1/20        # Natural death: ~1/age_max_bacteria per step

diffusion_speed_microbe = param_set_use$diffusion_speed_microbe
decay_rate_microbe      = param_set_use$decay_rate_microbe

# Macrophage pool diffusion (cells move slower than signals)
# diffusion_speed_macro = 0.05    # Very slow - cells don't diffuse much
diffusion_speed_macro = param_set_use$diffusion_speed_macro
decay_rate_macro = 0.0 # No natural decay - conserved pool

# Treg pool diffusion
diffusion_speed_treg = diffusion_speed_microbe     # Same as macrophages
decay_rate_treg = 0.0 # No natural decay - conserved pool

# Chemotactic sensitivity (Keller-Segel chi parameter)
# Macrophages chemotax toward danger signals (DAMPs + PAMPs)
# Tregs chemotax toward DAMPs
# Set chi = 0 to recover isotropic (homogeneous) diffusion
# 
# If your attractant is 0–1 and changes smoothly, maybe neighbor jumps are ~0.05–0.2
# χ_max safe ~ 2.5

chi_macro = 5
chi_treg  = 5

# Max density values (normalized to 1.0)
max_density_microbe = 1.0
max_density_macro = 1.0
max_density_treg = 1.0

active_age_limit = 1 # OVERWRITE FOR DIFFUSION MODEL!

density_M0_0 = 1
density_treg_0 = 1

extinction_limit = 1e-6 # 1e-3 is too big I think

rate_of_activation   = 0.05
rate_of_deactivation = 0.05


n_chemotaxis_smooth = 2 # 0 is no smoothing for the chemotaxis fields!

