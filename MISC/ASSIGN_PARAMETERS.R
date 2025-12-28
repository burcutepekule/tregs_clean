# ============================================================================
# ASSIGN PARAMETERS FROM CSV (WITH MACROPHAGE SPECIFICITY)
# ============================================================================
# Thresholds
th_ROS_microbe = param_set_use$th_ROS_microbe
th_ROS_epith_injury = param_set_use$th_ROS_epith_injury
epith_recovery_chance = param_set_use$epith_recovery_chance

# Diffusion speeds
diffusion_speed_DAMPs = param_set_use$diffusion_speed_DAMPs
diffusion_speed_SAMPs = param_set_use$diffusion_speed_SAMPs
diffusion_speed_PAMPs = param_set_use$diffusion_speed_PAMPs
diffusion_speed_ROS   = param_set_use$diffusion_speed_ROS

# Signal production - OLD
# if(ros_level==2){
#   add_ROS = 1
# }else{
#   add_ROS = ros_level*param_set_use$add_ROS # to see whether infection resolves without ROS 
# }

add_ROS = ros_level*param_set_use$add_ROS # ros_level=0, control

add_DAMPs = param_set_use$add_DAMPs
add_SAMPs = param_set_use$add_SAMPs
add_PAMPs = param_set_use$add_PAMPs

# Decay rates
ros_decay = param_set_use$ros_decay
DAMPs_decay = param_set_use$DAMPs_decay
SAMPs_decay = param_set_use$SAMPs_decay
PAMPs_decay = param_set_use$PAMPs_decay

# Activation thresholds
activation_threshold_danger = param_set_use$activation_threshold_danger
activation_threshold_SAMPs  = param_set_use$activation_threshold_SAMPs

# Engulfment activities
activity_engulf_M0_baseline = param_set_use$activity_engulf_M0_baseline
activity_engulf_M1_baseline = param_set_use$activity_engulf_M1_baseline
activity_engulf_M2_baseline = param_set_use$activity_engulf_M2_baseline

# ROS production activities
activity_ROS_M0_baseline = 0
activity_ROS_M2_baseline = 0
activity_ROS_M1_baseline = param_set_use$activity_ROS_M1_baseline

# Leak rates
rate_leak_commensal_injury   = param_set_use$rate_leak_commensal_injury
rate_leak_commensal_baseline = param_set_use$rate_leak_commensal_baseline
rate_leak_pathogen_injury = ifelse(sterile == 1, 0.0, param_set_use$rate_leak_pathogen_injury)
rate_leak_pathogen_injury = pat_level*rate_leak_pathogen_injury #1, 2, 3

# Phagocyte parameters
active_age_limit = as.integer(param_set_use$active_age_limit)
cc_phagocyte     = as.integer(param_set_use$cc_phagocyte)
# digestion_time   = 1 # by default digestion happens every time point

# Treg parameters
treg_vicinity_effect = 1
treg_discrimination_efficiency = param_set_use$treg_discrimination_efficiency
# allow_tregs_to_suppress_cognate = FALSE

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

