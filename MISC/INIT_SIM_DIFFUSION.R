# ============================================================================
# INITIALIZATION FOR DIFFUSION-BASED MODEL
# ============================================================================
# This file initializes all density grids and tracking matrices for the
# diffusion-based simulation model.
# ============================================================================

# ============================================================================
# INITIALIZE SIGNAL GRIDS (same as ABM)
# ============================================================================
DAMPs = matrix(0, grid_size, grid_size)
SAMPs = matrix(0, grid_size, grid_size)
PAMPs = matrix(0, grid_size, grid_size)
ROS   = matrix(0, grid_size, grid_size)

# ============================================================================
# INITIALIZE EPITHELIUM (same as ABM)
# ============================================================================
epithelium = data.frame(
  x = seq(1, grid_size, 1),
  y = rep(0, grid_size),
  level_injury = 0,
  id = seq(1, grid_size)
)
epithelium[injury_site, ]$level_injury = initial_injury_level

# ============================================================================
# INITIALIZE MICROBE DENSITY GRIDS
# ============================================================================
# Pathogen density: initial concentration at injury sites
density_pathogen = matrix(0, grid_size, grid_size)
if (n_pathogens_lp > 0) {
  # Distribute initial pathogen density across injury sites at epithelium (row 1)
  # Normalize to density (n_pathogens / area -> density per cell)
  initial_pathogen_density = n_pathogens_lp / length(injury_site)
  initial_pathogen_density = min(initial_pathogen_density, max_density_microbe)
  density_pathogen[1, injury_site] = initial_pathogen_density
}

# Commensal density: initial concentration spread across grid
density_commensal = matrix(0, grid_size, grid_size)
if (n_commensals_lp > 0) {
  # Distribute initial commensal density uniformly
  initial_commensal_density = n_commensals_lp / (grid_size * grid_size)
  initial_commensal_density = min(initial_commensal_density, max_density_microbe)
  density_commensal = matrix(initial_commensal_density, grid_size, grid_size)
}

# ============================================================================
# INITIALIZE MACROPHAGE DENSITY GRID
# ============================================================================
# Total macrophage pool (M0 + M1 + M2) - conserved quantity
# Initialize uniformly across the grid, normalized to density
density_macro = matrix(n_phagocytes / (grid_size * grid_size), grid_size, grid_size)
density_macro = pmin(density_macro, max_density_macro)

# M1 and M2 phenotype densities (initially zero - all are M0)
density_M1 = matrix(0, grid_size, grid_size)
density_M2 = matrix(0, grid_size, grid_size)

# ============================================================================
# INITIALIZE TREG DENSITY GRID
# ============================================================================
# Total Treg pool - conserved quantity
density_treg = matrix(n_tregs / (grid_size * grid_size), grid_size, grid_size)
density_treg = pmin(density_treg, max_density_treg)

# Active Treg density (initially zero - all are resting)
density_Treg_active = matrix(0, grid_size, grid_size)

# ============================================================================
# INITIALIZE DEATH COUNTERS
# ============================================================================
pathogens_killed_by_ROS = 0
commensals_killed_by_ROS = 0
# Note: No "killed by Mac" in diffusion model (no engulfment)

# ============================================================================
# INITIALIZE LONGITUDINAL TRACKING MATRICES
# ============================================================================
# Match the output format of the ABM for compatibility
epithelium_longitudinal  = matrix(0, nrow = t_max, ncol = 1)  # epithelial_score
macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3)  # M0, M1, M2
microbes_longitudinal    = matrix(0, nrow = t_max, ncol = 2)  # commensal, pathogen
tregs_longitudinal       = matrix(0, nrow = t_max, ncol = 2)  # resting, active
microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2*4)  # C_ROS,C_M0,C_M1,C_M2,P_ROS,P_M0,P_M1,P_M2

# Optional: grids for spatial tracking (if needed)
injury_pathogen_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
injury_ros_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
pathogen_epithelium_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)
ros_epithelium_longitudinal = matrix(0, nrow = t_max, ncol = grid_size)

# Initialize injury site tracker
injury_site_updated = injury_site
