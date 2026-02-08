# ============================================================================
# MAIN SIMULATION LOOP-DIFFUSION-BASED MODEL
# ============================================================================
# This replaces the agent-based model with a density-field approach where
# microbes and lymphocytes are represented as continuous density grids
# that evolve via diffusion, decay, and reaction terms.
#
# Key simplifications:
#-No individual agent tracking (no bacteria registry, no agent IDs)
#-All entities represented as density fields on the same 25x25 grid
#-Phenotype updates happen every active_age_limit steps
#-Much faster execution while preserving spatial dynamics
# ============================================================================

density_M0   = matrix(1, nrow = dim(density_M0)[1], ncol = dim(density_M0)[1])
density_treg = matrix(1, nrow = dim(density_treg)[1], ncol = dim(density_treg)[1])
# ============================================================================
# 0. EXTINCTION
# ============================================================================
epithelium[epithelium<extinction_limit]=0
density_pathogen[density_pathogen<extinction_limit] =0
density_commensal[density_commensal<extinction_limit]=0
DAMPs[DAMPs<extinction_limit] =0
SAMPs[SAMPs<extinction_limit]=0
PAMPs[PAMPs<extinction_limit]=0
ROS[ROS<extinction_limit]=0

density_M0[density_M0<extinction_limit]=0
density_M1[density_M1<extinction_limit]=0
density_M2[density_M2<extinction_limit]=0
density_treg[density_treg<extinction_limit]=0
density_treg_active[density_treg_active<extinction_limit]=0

# ============================================================================
# 1. MICROBE SOURCE TERMS (at epithelium, y=1)
# ============================================================================


# Pathogens leak through injured epithelium
pathogen_source       = epithelium$level_injury*rate_leak_pathogen_injury*0.01 # this is because I don't wanna change the rate_leak_pathogen_injury in ASSIGN_PARAMETERS.R
density_pathogen[1, ] = density_pathogen[1, ]+pathogen_source
density_pathogen[1, ] = pmin(density_pathogen[1, ], max_density_microbe)

# Commensals: baseline leak+injury-enhanced leak
commensal_source_baseline = rep(rate_leak_commensal_baseline*0.01)
commensal_source_injury   = epithelium$level_injury*rate_leak_commensal_injury*0.01
density_commensal[1, ] = density_commensal[1, ]+commensal_source_baseline+commensal_source_injury
density_commensal[1, ] = pmin(density_commensal[1, ], max_density_microbe)

# ============================================================================
# 2. DIFFUSE ALL DENSITY GRIDS
# ============================================================================
# Microbes
density_pathogen = diffuse_matrix_cpp(density_pathogen, diffusion_speed_microbe,
                                      max_density_microbe, reflect_top = FALSE)
density_commensal = diffuse_matrix_cpp(density_commensal, diffusion_speed_microbe,
                                       max_density_microbe, reflect_top = FALSE)

# Signals (same as before)
DAMPs = diffuse_matrix_cpp(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs, reflect_top = FALSE)
SAMPs = diffuse_matrix_cpp(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs, reflect_top = FALSE)
PAMPs = diffuse_matrix_cpp(PAMPs, diffusion_speed_PAMPs, max_cell_value_PAMPs, reflect_top = FALSE)
ROS   = diffuse_matrix_cpp(ROS, diffusion_speed_ROS, max_cell_value_ROS, reflect_top = FALSE)


# Macrophage pools: chemotax toward danger signals (DAMPs + PAMPs)
danger_signal_grid = DAMPs + PAMPs

# Smooth the chemotactic field so the gradient extends further into tissue.
# The raw DAMPs+PAMPs field decays too steeply from the epithelium
# (exp(-x*sqrt(decay/D)) ≈ near-zero by row 3-4), concentrating the
# chemotactic gradient in 1-2 rows and creating banding artifacts.
# Extra diffusion passes spread the gradient without changing actual signals.
chemotaxis_field_macro = danger_signal_grid
if(n_chemotaxis_smooth>0){
  for (s in seq_len(n_chemotaxis_smooth)) {
    chemotaxis_field_macro = diffuse_matrix_cpp(chemotaxis_field_macro, 0.12,
                                                max(chemotaxis_field_macro) + 1e-10,
                                                reflect_top = TRUE)
  }
}


## NO NEED TO DIFFUSE density_M0, ALWAYS 1 EVERYWHERE
# density_M0   = diffuse_matrix_biased_cpp(density_M0, chemotaxis_field_macro,
#                                          diffusion_speed_macro, chi_macro,
#                                          max_density_macro, reflect_top = TRUE)
density_M1   = diffuse_matrix_biased_cpp(density_M1, chemotaxis_field_macro,
                                         diffusion_speed_macro, chi_macro,
                                         max_density_macro, reflect_top = TRUE)
density_M2   = diffuse_matrix_biased_cpp(density_M2, chemotaxis_field_macro,
                                         diffusion_speed_macro, chi_macro,
                                         max_density_macro, reflect_top = TRUE)

if(randomize_tregs==1){
  
  density_treg_active = diffuse_matrix_cpp(density_treg_active, diffusion_speed_treg, max_density_treg, reflect_top = TRUE)
  
}else{
  # Treg pools: chemotax toward DAMPs (also smoothed)
  chemotaxis_field_treg = DAMPs
  if(n_chemotaxis_smooth>0){
    for (s in seq_len(n_chemotaxis_smooth)) {
      chemotaxis_field_treg = diffuse_matrix_cpp(chemotaxis_field_treg, 0.12,
                                                 max(chemotaxis_field_treg) + 1e-10,
                                                 reflect_top = TRUE)}
  }
  
  # density_treg = diffuse_matrix_biased_cpp(density_treg, DAMPs,
  #                                          diffusion_speed_treg, chi_treg,
  #                                          max_density_treg, reflect_top = TRUE)
  density_treg_active = diffuse_matrix_biased_cpp(density_treg_active, DAMPs,
                                                  diffusion_speed_treg, chi_treg,
                                                  max_density_treg, reflect_top = TRUE)
}


# ============================================================================
# 3. DECAY ALL GRIDS
# ============================================================================
# Microbe decay (natural death)
density_pathogen  = density_pathogen*(1-decay_rate_microbe)
density_commensal = density_commensal*(1-decay_rate_microbe)

# # Macrophage and Treg pools (no decay if conserved)
# density_macro = density_macro*(1-decay_rate_macro)
# density_treg  = density_treg*(1-decay_rate_treg)

# Signal decay
DAMPs = DAMPs*(1-DAMPs_decay)
SAMPs = SAMPs*(1-SAMPs_decay)
PAMPs = PAMPs*(1-PAMPs_decay)
ROS   = ROS*(1-ros_decay)

# ============================================================================
# 4. SIGNAL PRODUCTION
# ============================================================================
# DAMPs from epithelial injury (at y=1)
DAMPs[1, ] = DAMPs[1, ]+epithelium$level_injury*add_DAMPs

# DAMPs from microbes at epithelium (both pathogens and commensals cause tissue stress)
microbe_density_at_epith = density_pathogen[1, ]+density_commensal[1, ]
DAMPs[1, ]               = DAMPs[1, ]+microbe_density_at_epith*add_DAMPs
DAMPs                    = pmin(DAMPs, max_cell_value_DAMPs)

# PAMPs from pathogen density (throughout the grid)
PAMPs = PAMPs+density_pathogen*add_PAMPs
PAMPs = pmin(PAMPs, max_cell_value_PAMPs)

# ROS from M1 macrophage density
ROS = ROS+density_M1*activity_ROS_M1_baseline*add_ROS
ROS = pmin(ROS, max_cell_value_ROS)

# SAMPs from active Treg density
SAMPs = SAMPs+density_treg_active*add_SAMPs*allow_tregs
SAMPs = pmin(SAMPs, max_cell_value_SAMPs)

# ============================================================================
# 5a. KILL MICROBES WITH ROS
# ============================================================================
# Microbes die where ROS exceeds threshold
# Using smooth killing: kill rate proportional to how much ROS exceeds threshold
ros_kill_factor = pmax(0, (ROS-th_ROS_microbe) / (1-th_ROS_microbe))

# ============================================================================
# 5b. KILL MICROBES WITH ENGULFMENT
# ============================================================================

eng_kill_factor_M1 = density_M1*activity_engulf_M1_baseline 
eng_kill_factor_M2 = density_M2*activity_engulf_M2_baseline
eng_kill_factor = eng_kill_factor_M1 + eng_kill_factor_M2

total_kill_factor = pmin(ros_kill_factor+eng_kill_factor, 1)  # Cap at 1

# (when calculated like below sum might be more than 1 but ratios will be used so should be fine)
# Track killed microbes for longitudinal data
pathogens_killed_by_ROS_this_step = sum(density_pathogen*ros_kill_factor)
commensals_killed_by_ROS_this_step = sum(density_commensal*ros_kill_factor)

# Track engulfed microbes for ANTIGEN PRESENTATION
pathogens_engf_M1  = density_pathogen*eng_kill_factor_M1
pathogens_engf_M2  = density_pathogen*eng_kill_factor_M2
commensals_engf_M1 = density_commensal*eng_kill_factor_M1
commensals_engf_M2 = density_commensal*eng_kill_factor_M2

# Apply killing
density_pathogen  = density_pathogen*(1-total_kill_factor)
density_commensal = density_commensal*(1-total_kill_factor)

# Update cumulative counters
pathogens_killed_by_ROS = pathogens_killed_by_ROS+pathogens_killed_by_ROS_this_step
commensals_killed_by_ROS = commensals_killed_by_ROS+commensals_killed_by_ROS_this_step

# ============================================================================
# 6. UPDATE EPITHELIAL INJURY
# ============================================================================
# Injury from pathogens at epithelium
pathogen_injury = density_pathogen[1, ]*c_in_log  # Scale factor for injury
epithelium$level_injury = epithelium$level_injury+pathogen_injury

# Injury from ROS (ROS above threshold damages epithelium)
ros_at_epith = ROS[1, ]
ros_injury = as.numeric(ros_at_epith > th_ROS_epith_injury)
epithelium$level_injury = epithelium$level_injury+ros_injury

# Cap injury at maximum
epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)

# Stochastic recovery
recovery_mask = runif(grid_size) < epith_recovery_chance
epithelium$level_injury = epithelium$level_injury-recovery_mask*(epithelium$level_injury > 0)
epithelium$level_injury = pmax(epithelium$level_injury, 0)

# Update injury site
injury_site_updated = which(epithelium$level_injury > 0)

# ============================================================================
# 7. PHENOTYPE UPDATES (every active_age_limit steps)
# ============================================================================
# if (t %% active_age_limit == 0) {

# Compute danger signal grid
danger_signal_grid = DAMPs+PAMPs

DANGER_diff=matrix(pmax(0, danger_signal_grid - activation_threshold_danger*1e-3)/ activation_threshold_danger*1e-3, 
                   nrow = nrow(danger_signal_grid), 
                   ncol = ncol(danger_signal_grid))

SAMPS_diff=matrix(pmax(0, SAMPs - activation_threshold_SAMPs*1e-3) / activation_threshold_SAMPs*1e-3, 
                  nrow = nrow(SAMPs), 
                  ncol = ncol(SAMPs))

# Total activation signal
total_diff = DANGER_diff+SAMPS_diff

# Soft split: fraction that becomes M1 vs M2
# Where total_diff > 0, compute proportional split
# Where total_diff == 0, no activation (M0)
frac_M1_add = safe_divide(DANGER_diff, total_diff)
frac_M2_add = safe_divide(SAMPS_diff, total_diff)

# # when both signals are off
# frac_M_remove = floor(1-frac_M1_add) + floor(1-frac_M2_add)
# frac_M_remove = matrix(ifelse(frac_M_remove==2,1,0), 
#                        nrow = nrow(frac_M_remove), 
#                        ncol = ncol(frac_M_remove))
# frac_M_remove+frac_M1_add+frac_M2_add should be 1!

# when both signals are off
frac_M1_remove = floor(1-frac_M1_add)
frac_M2_remove = floor(1-frac_M2_add)

# avoid depleting all at the same time when frac=1
d_frac_M1 = density_M0*rate_of_activation*frac_M1_add - density_M1*rate_of_deactivation*frac_M1_remove - density_M1*rate_of_activation*frac_M2_add
d_frac_M2 = density_M1*rate_of_activation*frac_M2_add - density_M2*rate_of_deactivation*frac_M2_remove

# Apply fractions to macrophage pool
density_M1 = density_M1+d_frac_M1
density_M2 = density_M2+d_frac_M2

# Cap densities
density_M1 = pmin(density_M1, max_density_macro)
density_M2 = pmin(density_M2, max_density_macro)

# recalculate remaining M0
# density_M0 = density_M0 - (density_M1 + density_M2) # here all of them are turning so 

# ========================================================================
# TREG ACTIVATION
# ========================================================================
# Tregs activate where macrophages (M1 or M2) are present
# Activation probability depends on commensal/pathogen ratio

pathogen_density_M1  = safe_divide(pathogens_engf_M1,(pathogens_engf_M1+commensals_engf_M1))
commensal_density_M1 = safe_divide(commensals_engf_M1,(pathogens_engf_M1+commensals_engf_M1))
ag_presentation_commensal_M1 = safe_divide(commensal_density_M1,(pathogen_density_M1+commensal_density_M1))

pathogen_density_M2  = safe_divide(pathogens_engf_M2,(pathogens_engf_M2+commensals_engf_M2))
commensal_density_M2 = safe_divide(commensals_engf_M2,(pathogens_engf_M2+commensals_engf_M2))
ag_presentation_commensal_M2 = safe_divide(commensal_density_M2,(pathogen_density_M2+commensal_density_M2))

treg_activation_by_M1_prob = sample_with_efficiency_beta(ag_presentation_commensal_M1, treg_discrimination_efficiency)

# treg_activation_by_M1_prob = ag_presentation_commensal_M1*treg_discrimination_efficiency +
#   (1-ag_presentation_commensal_M1)*(1-treg_discrimination_efficiency)

# treg_activation_by_M2_prob = ag_presentation_commensal_M2*treg_discrimination_efficiency +
#   (1-ag_presentation_commensal_M2)*(1-treg_discrimination_efficiency)

treg_activation_by_M2_prob = 0

treg_activation_prob = treg_activation_by_M1_prob + treg_activation_by_M2_prob
treg_activation_prob = matrix(pmin(1, treg_activation_prob), 
                              nrow = nrow(treg_activation_prob), 
                              ncol = ncol(treg_activation_prob))

d_tregs_active =  rate_of_activation*density_treg*treg_activation_prob - density_treg_active*rate_of_deactivation
# Update active Treg density
density_treg_active = d_tregs_active*allow_tregs # avoid depleting all when treg_activation_prob=1

density_treg_active = pmin(density_treg_active, max_density_treg)
# }

# ============================================================================
# 0. EXTINCTION
# ============================================================================
epithelium[epithelium<extinction_limit]=0
density_pathogen[density_pathogen<extinction_limit] =0
density_commensal[density_commensal<extinction_limit]=0
DAMPs[DAMPs<extinction_limit] =0
SAMPs[SAMPs<extinction_limit]=0
PAMPs[PAMPs<extinction_limit]=0
ROS[ROS<extinction_limit]=0

density_M0[density_M0<extinction_limit]=0
density_M1[density_M1<extinction_limit]=0
density_M2[density_M2<extinction_limit]=0
density_treg[density_treg<extinction_limit]=0
density_treg_active[density_treg_active<extinction_limit]=0

# ============================================================================
# 8. SAVE LONGITUDINAL DATA
# ============================================================================
# Epithelial score (higher = more injured)
epithelial_score = sum(epithelium$level_injury)
epithelium_longitudinal[t, 1] = epithelial_score

# Macrophage densities (sum across grid, then convert to "counts")
# M0 = total macro pool-M1-M2
total_M1    = sum(density_M1)
total_M2    = sum(density_M2)
total_M0    = sum(density_M0)

macrophages_longitudinal[t, ] = c(total_M0, total_M1, total_M2)

# Microbe densities (sum across grid, scaled)
total_commensal = sum(density_commensal)
total_pathogen = sum(density_pathogen)
microbes_longitudinal[t, ] = c(total_commensal, total_pathogen)

# Treg densities
total_treg = sum(density_treg)
total_treg_active = sum(density_treg_active)
total_treg_resting = max(0, total_treg-total_treg_active)
# tregs_longitudinal[t, ] = c(total_treg_resting, total_treg_active)*grid_size^2 / max_density_treg
tregs_longitudinal[t, ] = c(total_treg_resting, total_treg_active)

# Cumulative death tracking (scaled for comparability)
# Note: In diffusion model, "killed by Mac" doesn't apply-only ROS killing
microbes_cumdeath_longitudinal[t, ] = c(
  commensals_killed_by_ROS, 0, 0, 0,  # C_ROS, C_M0, C_M1, C_M2
  pathogens_killed_by_ROS, 0, 0, 0    # P_ROS, P_M0, P_M1, P_M2
)

# Add this at the end of your simulation loop (replacing the image() call)

if(plot_grid_t==1){
  # Function to plot matrix in "natural" orientation
  plot_grid = function(mat, title, zlim = c(0, 0.5)) {
    image(
      t(mat)[, nrow(mat):1],
      main = title,
      zlim = zlim,
      col = hcl.colors(50, "viridis"),
      axes = FALSE
    )
  }
  
  # Create faceted plot
  png(filename = paste0(dir_name_data,"/grids_",t,".png"), width = 1200, height = 800)
  par(mfrow = c(3, 4), mar = c(2, 2, 3, 1))
  
  plot_grid(density_pathogen, "Pathogen", zlim = c(0, 1.01*max(density_pathogen)))
  plot_grid(density_commensal, "Commensal", zlim = c(0, 1.01*max(density_commensal)))
  plot_grid(density_M0, "M0", zlim = c(0, 1.01*max(density_M0)))
  plot_grid(density_M1, "M1", zlim = c(0, 1.01*max(density_M1)))
  plot_grid(density_M2, "M2", zlim = c(0, 1.01*max(density_M2)))
  plot_grid(density_treg, "Treg (resting)", zlim = c(0, 1.01*max(density_treg)))
  plot_grid(density_treg_active, "Treg (active)", zlim = c(0, 1.01*max(density_treg_active)))
  
  plot_grid(DAMPs, "DAMPs", zlim = c(0, 1.01*max(DAMPs)))
  plot_grid(PAMPs, "PAMPs", zlim = c(0, 1.01*max(PAMPs)))
  plot_grid(SAMPs, "SAMPs", zlim = c(0, 1.01*max(SAMPs)))
  plot_grid(ROS, "ROS", zlim = c(0, 1.01*max(ROS)))
  
  # Epithelial injury as a 1D barplot
  barplot(epithelium$level_injury, main = "Epithelial injury", 
          col = "darkblue", border = NA, ylim = c(0, max_level_injury))
  
  
  # Add overall title
  mtext(sprintf("Time: %d", t), outer = TRUE, line = -1.5, cex = 1.5, font = 2)
  
  dev.off()
}

