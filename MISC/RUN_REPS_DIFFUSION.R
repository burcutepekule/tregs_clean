# ============================================================================
# DIFFUSION MODEL - REPLICATE RUNNER
# ============================================================================
# This is the equivalent of RUN_REPS_CPP_ABM_PAMPS.R but for the
# diffusion-based model. Runs num_reps replicates.
# ============================================================================

for (reps_in in 0:(num_reps-1)){

  # ========================================================================
  # INITIALIZE SIGNAL GRIDS
  # ========================================================================
  DAMPs = matrix(0, grid_size, grid_size)
  SAMPs = matrix(0, grid_size, grid_size)
  PAMPs = matrix(0, grid_size, grid_size)
  ROS   = matrix(0, grid_size, grid_size)

  # ========================================================================
  # INITIALIZE EPITHELIUM
  # ========================================================================
  epithelium = data.frame(
    x = seq(1, grid_size, 1),
    y = rep(0, grid_size),
    level_injury = 0,
    id = seq(1, grid_size)
  )
  epithelium[injury_site, ]$level_injury = initial_injury_level

  # ========================================================================
  # INITIALIZE MICROBE DENSITY GRIDS
  # ========================================================================
  # Pathogen density: initial concentration at injury sites at epithelium (row 1)
  density_pathogen = matrix(0, grid_size, grid_size)
  # if (n_pathogens_lp > 0) {
  #   initial_pathogen_density = n_pathogens_lp / length(injury_site)
  #   initial_pathogen_density = min(initial_pathogen_density, max_density_microbe)
  #   density_pathogen[1, injury_site] = initial_pathogen_density
  # }

  # Commensal density: initial concentration spread across grid
  density_commensal = matrix(0, grid_size, grid_size)
  # if (n_commensals_lp > 0) {
  #   initial_commensal_density = n_commensals_lp / (grid_size * grid_size)
  #   initial_commensal_density = min(initial_commensal_density, max_density_microbe)
  #   density_commensal = matrix(initial_commensal_density, grid_size, grid_size)
  # }

  # ========================================================================
  # INITIALIZE MACROPHAGE DENSITY GRID
  # ========================================================================
  # Total macrophage pool (M0 + M1 + M2) - conserved quantity
  density_macro = matrix(n_phagocytes / (grid_size * grid_size), grid_size, grid_size)
  density_macro = pmin(density_macro, max_density_macro)

  # M1 and M2 phenotype densities (initially zero - all are M0)
  density_M1 = matrix(0, grid_size, grid_size)
  density_M2 = matrix(0, grid_size, grid_size)

  # ========================================================================
  # INITIALIZE TREG DENSITY GRID
  # ========================================================================
  density_treg = matrix(n_tregs / (grid_size * grid_size), grid_size, grid_size)
  density_treg = pmin(density_treg, max_density_treg)

  # Active Treg density (initially zero - all are resting)
  density_Treg_active = matrix(0, grid_size, grid_size)

  # ========================================================================
  # INITIALIZE DEATH COUNTERS
  # ========================================================================
  pathogens_killed_by_ROS = 0
  commensals_killed_by_ROS = 0

  # ========================================================================
  # INITIALIZE LONGITUDINAL TRACKING MATRICES
  # ========================================================================
  epithelium_longitudinal  = matrix(0, nrow = t_max, ncol = 1)
  macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3)
  microbes_longitudinal    = matrix(0, nrow = t_max, ncol = 2)
  tregs_longitudinal       = matrix(0, nrow = t_max, ncol = 2)
  microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2*4)
  injury_pathogen_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
  injury_ros_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
  pathogen_epithelium_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
  ros_epithelium_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)

  # Initialize injury site tracker
  injury_site_updated = injury_site

  # ========================================================================
  # MAIN SIMULATION LOOP
  # ========================================================================
  source('./MISC/MAIN_SIM_LOOP_OVER_T_DIFFUSION.R')

  # ========================================================================
  # SAVE LONGITUDINAL DATA
  # ========================================================================
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
  longitudinal_df$rep_id = reps_in

  # Compute steady state indices
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
                  rep_id,
                  epithelial_score,
                  time_ss_e,
                  time_ss_p,
                  everything())

  longitudinal_df_keep = rbind(longitudinal_df_keep, longitudinal_df)
}
