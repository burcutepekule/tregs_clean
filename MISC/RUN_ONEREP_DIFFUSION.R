# ============================================================================
# DIFFUSION MODEL - SINGLE REPLICATE RUNNER
# ============================================================================
# This is the equivalent of RUN_ONEREP_CPP_ABM_PAMPS.R but for the
# diffusion-based model.
# ============================================================================

# Initialize all density grids and tracking matrices
source('./MISC/INIT_SIM_DIFFUSION.R')

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
source('./MISC/MAIN_SIM_LOOP_OVER_T_DIFFUSION.R')

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
