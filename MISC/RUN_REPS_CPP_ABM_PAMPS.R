# ============================================================================
# C++ ACCELERATED VERSION WITH PAMPs SUPPORT
# ============================================================================

for (reps_in in 0:(num_reps-1)){
  
  # Initialize fields
  DAMPs = matrix(0, grid_size, grid_size)
  SAMPs = matrix(0, grid_size, grid_size)
  PAMPs = matrix(0, grid_size, grid_size)  # Toxins released by pathogens that diffuse through the environment
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
  phagocyte_phenotype = rep(0, n_phagocytes)  # 0=M0, 1=M1, 2=M2
  phagocyte_activity_ROS = rep(activity_ROS_M0_baseline, n_phagocytes)
  phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
  phagocyte_active_age = rep(0, n_phagocytes)
  phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)  # Memory: +1=commensal, -1=pathogen, 0=empty
  
  # Initialize tregs
  treg_x = sample(1:grid_size, n_tregs, TRUE)
  treg_y = sample(1:grid_size, n_tregs, TRUE)
  treg_active_age = rep(0, n_tregs)
  treg_phenotype = rep(0, n_tregs)  # 0=resting, 1=activated
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
  macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3) #onelevel
  microbes_longitudinal    = matrix(0, nrow = t_max, ncol = 2)
  tregs_longitudinal       = matrix(0, nrow = t_max, ncol = 2)
  microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2*4)
  
  # ============================================================================
  # MAIN SIMULATION LOOP
  # ============================================================================
  source('./MISC/MAIN_SIM_LOOP_OVER_T.R')
  
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
  
  
  longitudinal_df = longitudinal_df %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                                        5*epithelial_inj_1+
                                                        4*epithelial_inj_2+
                                                        3*epithelial_inj_3+
                                                        2*epithelial_inj_4+
                                                        1*epithelial_inj_5)
  
  # moved from one time_ss based only on epithelial_score to two time_ss for both signals
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
