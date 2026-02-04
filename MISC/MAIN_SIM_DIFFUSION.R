# ============================================================================
# MAIN SIMULATION LOOP - DIFFUSION-BASED MODEL
# ============================================================================
# This replaces the agent-based model with a density-field approach where
# microbes and lymphocytes are represented as continuous density grids
# that evolve via diffusion, decay, and reaction terms.
#
# Key simplifications:
# - No individual agent tracking (no bacteria registry, no agent IDs)
# - All entities represented as density fields on the same 25x25 grid
# - Phenotype updates happen every active_age_limit steps
# - Much faster execution while preserving spatial dynamics
# ============================================================================

# ============================================================================
# 1. MICROBE SOURCE TERMS (at epithelium, y=1)
# ============================================================================
# Pathogens leak through injured epithelium
pathogen_source       = epithelium$level_injury*rate_leak_pathogen_injury*0.01 # this is because I don't wanna change the rate_leak_pathogen_injury in ASSIGN_PARAMETERS.R
density_pathogen[1, ] = density_pathogen[1, ] + pathogen_source
density_pathogen[1, ] = pmin(density_pathogen[1, ], max_density_microbe)

# Commensals: baseline leak + injury-enhanced leak
commensal_source_baseline = rep(rate_leak_commensal_baseline*0.01)
commensal_source_injury   = epithelium$level_injury*rate_leak_commensal_injury*0.01
density_commensal[1, ] = density_commensal[1, ] + commensal_source_baseline + commensal_source_injury
density_commensal[1, ] = pmin(density_commensal[1, ], max_density_microbe)

# ============================================================================
# 2. DIFFUSE ALL DENSITY GRIDS
# ============================================================================
# Microbes
density_pathogen = diffuse_matrix_cpp(density_pathogen, diffusion_speed_microbe,
                                       max_density_microbe, reflect_top = FALSE)
density_commensal = diffuse_matrix_cpp(density_commensal, diffusion_speed_microbe,
                                        max_density_microbe, reflect_top = FALSE)

# Macrophage and Treg pools
density_macro = diffuse_matrix_cpp(density_macro, diffusion_speed_macro,
                                    max_density_macro, reflect_top = FALSE)
density_treg = diffuse_matrix_cpp(density_treg, diffusion_speed_treg,
                                   max_density_treg, reflect_top = FALSE)

# Signals (same as before)
DAMPs = diffuse_matrix_cpp(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs, reflect_top = FALSE)
SAMPs = diffuse_matrix_cpp(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs, reflect_top = FALSE)
PAMPs = diffuse_matrix_cpp(PAMPs, diffusion_speed_PAMPs, max_cell_value_PAMPs, reflect_top = FALSE)
ROS   = diffuse_matrix_cpp(ROS, diffusion_speed_ROS, max_cell_value_ROS, reflect_top = FALSE)

# ============================================================================
# 3. DECAY ALL GRIDS
# ============================================================================
# Microbe decay (natural death)
density_pathogen  = density_pathogen*(1-decay_rate_microbe)
density_commensal = density_commensal*(1-decay_rate_microbe)

# # Macrophage and Treg pools (no decay if conserved)
# density_macro = density_macro * (1-decay_rate_macro)
# density_treg  = density_treg * (1 - decay_rate_treg)

# Signal decay
DAMPs = DAMPs * (1-DAMPs_decay)
SAMPs = SAMPs * (1-SAMPs_decay)
PAMPs = PAMPs * (1-PAMPs_decay)
ROS   = ROS * (1-ros_decay)

# ============================================================================
# 4. SIGNAL PRODUCTION
# ============================================================================
# DAMPs from epithelial injury (at y=1)
DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury*add_DAMPs

# DAMPs from microbes at epithelium (both pathogens and commensals cause tissue stress)
microbe_density_at_epith = density_pathogen[1, ] + density_commensal[1, ]
DAMPs[1, ] = DAMPs[1, ] + microbe_density_at_epith*add_DAMPs
DAMPs = pmin(DAMPs, max_cell_value_DAMPs)

# PAMPs from pathogen density (throughout the grid)
PAMPs = PAMPs + density_pathogen*add_PAMPs
PAMPs = pmin(PAMPs, max_cell_value_PAMPs)

# ROS from M1 macrophage density
ROS = ROS + density_M1*activity_ROS_M1_baseline * add_ROS
ROS = pmin(ROS, max_cell_value_ROS)

# SAMPs from active Treg density
SAMPs = SAMPs + density_Treg_active * add_SAMPs * allow_tregs
SAMPs = pmin(SAMPs, max_cell_value_SAMPs)

# ============================================================================
# 5. KILL MICROBES WITH ROS
# ============================================================================
# Microbes die where ROS exceeds threshold
# Using smooth killing: kill rate proportional to how much ROS exceeds threshold
ros_kill_factor = pmax(0, (ROS - th_ROS_microbe) / (1 - th_ROS_microbe))
ros_kill_factor = pmin(ros_kill_factor, 1)  # Cap at 1

# Track killed microbes for longitudinal data
pathogens_killed_by_ROS_this_step = sum(density_pathogen * ros_kill_factor)
commensals_killed_by_ROS_this_step = sum(density_commensal * ros_kill_factor)

# Apply killing
density_pathogen = density_pathogen * (1 - ros_kill_factor)
density_commensal = density_commensal * (1 - ros_kill_factor)

# Update cumulative counters
pathogens_killed_by_ROS = pathogens_killed_by_ROS + pathogens_killed_by_ROS_this_step
commensals_killed_by_ROS = commensals_killed_by_ROS + commensals_killed_by_ROS_this_step

# ============================================================================
# 6. UPDATE EPITHELIAL INJURY
# ============================================================================
# Injury from pathogens at epithelium
pathogen_injury = density_pathogen[1, ]*c_in_log  # Scale factor for injury
epithelium$level_injury = epithelium$level_injury + pathogen_injury

# Injury from ROS (ROS above threshold damages epithelium)
ros_at_epith = ROS[1, ]
ros_injury = as.numeric(ros_at_epith > th_ROS_epith_injury)
epithelium$level_injury = epithelium$level_injury + ros_injury

# Cap injury at maximum
epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)

# Stochastic recovery
recovery_mask = runif(grid_size) < epith_recovery_chance
epithelium$level_injury = epithelium$level_injury - recovery_mask * (epithelium$level_injury > 0)
epithelium$level_injury = pmax(epithelium$level_injury, 0)

# Update injury site
injury_site_updated = which(epithelium$level_injury > 0)

# ============================================================================
# 7. PHENOTYPE UPDATES (every active_age_limit steps)
# ============================================================================
if (t %% active_age_limit == 0) {

  # Compute danger signal grid
  danger_signal_grid = DAMPs + PAMPs

  # Compute normalized differences from thresholds
  DANGER_diff = pmax(0, danger_signal_grid - activation_threshold_danger) / activation_threshold_danger
  SAMPS_diff = pmax(0, SAMPs - activation_threshold_SAMPs) / activation_threshold_SAMPs

  # Total activation signal
  total_diff = DANGER_diff + SAMPS_diff

  # Soft split: fraction that becomes M1 vs M2
  # Where total_diff > 0, compute proportional split
  # Where total_diff == 0, no activation (M0)
  frac_M1 = ifelse(total_diff > 0, DANGER_diff / total_diff, 0)
  frac_M2 = ifelse(total_diff > 0, SAMPS_diff / total_diff, 0)

  # Apply fractions to macrophage pool
  # Only activate where there's sufficient signal (at least one diff > 0)
  activation_mask = (total_diff > 0)

  density_M1 = density_macro * frac_M1 * activation_mask
  density_M2 = density_macro * frac_M2 * activation_mask

  # Cap densities
  density_M1 = pmin(density_M1, max_density_macro)
  density_M2 = pmin(density_M2, max_density_macro)

  # ========================================================================
  # TREG ACTIVATION
  # ========================================================================
  # Tregs activate where macrophages (M1 or M2) are present
  # Activation probability depends on commensal/pathogen ratio

  macro_presence = density_M1 + density_M2

  # Compute commensal ratio at each location
  total_microbe = density_commensal + density_pathogen
  commensal_ratio = ifelse(total_microbe > 0,
                           density_commensal / total_microbe,
                           0.5)  # Default to 0.5 if no microbes

  # Treg activation probability based on discrimination efficiency
  # P(activate) = P(commensal) * P(correctly identify) + P(pathogen) * P(incorrectly identify)
  # = ratio * efficiency + (1-ratio) * (1-efficiency)
  treg_activation_prob = commensal_ratio * treg_discrimination_efficiency +
                         (1 - commensal_ratio) * (1 - treg_discrimination_efficiency)

  # Scale by macrophage presence (need antigen presentation)
  treg_activation_prob = treg_activation_prob * pmin(macro_presence, 1)

  # Update active Treg density
  density_Treg_active = density_treg * treg_activation_prob * allow_tregs
  density_Treg_active = pmin(density_Treg_active, max_density_treg)
}

# ============================================================================
# 8. SAVE LONGITUDINAL DATA
# ============================================================================
# Epithelial score (higher = healthier)
# Score = sum of (max_injury - actual_injury) for each cell, normalized
epithelial_score = sum(max_level_injury - epithelium$level_injury)
epithelium_longitudinal[t, 1] = epithelial_score

# Macrophage densities (sum across grid, then convert to "counts")
# M0 = total macro pool - M1 - M2
total_M1 = sum(density_M1)
total_M2 = sum(density_M2)
total_macro = sum(density_macro)
total_M0 = max(0, total_macro - total_M1 - total_M2)

macrophages_longitudinal[t, ] = c(total_M0, total_M1, total_M2)/grid_size^2

# Microbe densities (sum across grid, scaled)
total_commensal = sum(density_commensal)
total_pathogen = sum(density_pathogen)
microbes_longitudinal[t, ] = c(total_commensal, total_pathogen)/grid_size^2

# Treg densities
total_treg = sum(density_treg)
total_treg_active = sum(density_Treg_active)
total_treg_resting = max(0, total_treg - total_treg_active)
# tregs_longitudinal[t, ] = c(total_treg_resting, total_treg_active) * grid_size^2 / max_density_treg
tregs_longitudinal[t, ] = c(total_treg_resting, total_treg_active)/grid_size^2

# Cumulative death tracking (scaled for comparability)
# Note: In diffusion model, "killed by Mac" doesn't apply - only ROS killing
microbes_cumdeath_longitudinal[t, ] = c(
  commensals_killed_by_ROS, 0, 0, 0,  # C_ROS, C_M0, C_M1, C_M2
  pathogens_killed_by_ROS, 0, 0, 0    # P_ROS, P_M0, P_M1, P_M2
)
