# ============================================================================
# C++ ACCELERATED VERSION WITH PAMPs SUPPORT + RANDOM SEARCH OPTIMIZATION
# ============================================================================

# Helper function to reset simulation state
reset_simulation_state = function() {
  DAMPs <<- matrix(0, grid_size, grid_size)
  SAMPs <<- matrix(0, grid_size, grid_size)
  PAMPs <<- matrix(0, grid_size, grid_size)
  ROS   <<- matrix(0, grid_size, grid_size)
  density_pathogen   <<- matrix(0, grid_size, grid_size)
  density_commensal  <<- matrix(0, grid_size, grid_size)
  density_M0  <<- matrix(0, grid_size, grid_size)
  density_M1  <<- matrix(0, grid_size, grid_size)
  density_M2  <<- matrix(0, grid_size, grid_size)
  density_treg  <<- matrix(0, grid_size, grid_size)
  density_treg_active  <<- matrix(0, grid_size, grid_size)
  
  epithelium$level_injury <<- 0
  epithelium[injury_site, ]$level_injury <<- 1
  
}

# Initialize fields
DAMPs = matrix(0, grid_size, grid_size)
SAMPs = matrix(0, grid_size, grid_size)
PAMPs = matrix(0, grid_size, grid_size)
ROS   = matrix(0, grid_size, grid_size)
density_pathogen = matrix(0, grid_size, grid_size)
density_commensal = matrix(0, grid_size, grid_size)
density_M0 = matrix(0, grid_size, grid_size)
density_M1 = matrix(0, grid_size, grid_size)
density_M2 = matrix(0, grid_size, grid_size)
density_treg = matrix(0, grid_size, grid_size)
density_treg_active = matrix(0, grid_size, grid_size)

# Initialize epithelium
epithelium = data.frame(
  x = seq(1, grid_size, 1),
  y = rep(0, grid_size),
  level_injury = 0,
  id = seq(1, grid_size)
)
epithelium[injury_site, ]$level_injury = 1

# Death counters
pathogens_killed_by_ROS = 0
pathogens_killed_by_Mac = rep(0, 3)
commensals_killed_by_ROS = 0
commensals_killed_by_Mac = rep(0, 3)

# Longitudinal tracking
epithelium_longitudinal  = matrix(0, nrow = t_max, ncol = 1)
macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3)
microbes_longitudinal    = matrix(0, nrow = t_max, ncol = 2)
tregs_longitudinal       = matrix(0, nrow = t_max, ncol = 2)
microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2*4)
pathogens_lumen_longitudinal   = matrix(0, nrow = t_max, ncol = 1)

# ============================================================================
# RANDOM SEARCH CONFIGURATION
# ============================================================================
param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

lower_bounds = param_bounds$lower
upper_bounds = param_bounds$upper

# Initialize success log file
success_log_file = paste0("./random_search_successes_", param_set_id_use, "_scenario_", scenario_ind, "_n2_", n2, ".txt")
# cat(paste0("# Random Search Success Log for ",param_set_id_use,"_pat_", pat_level, "_ros_", ros_level, "\n"),
#     file = success_log_file, append = FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("RANDOM SEARCH OPTIMIZATION STARTING\n"))
cat(sprintf("========================================\n"))

# ============================================================================
# MAIN RANDOM SEARCH LOOP
# ============================================================================
for (iter in 1:max_iterations) {
  
  # Sample random parameters
  current_theta = sapply(1:length(param_names), function(i) {
    runif(1, lower_bounds[i], upper_bounds[i])
  })
  current_theta = round(current_theta, 3)
  
  # Set parameter values
  diffusion_speed_SAMPs          = current_theta[1]
  add_SAMPs                      = current_theta[2]
  SAMPs_decay                    = current_theta[3]
  treg_discrimination_efficiency = current_theta[4]
  activation_threshold_SAMPs     = current_theta[5]
  
  cat(sprintf("\n>>> [ITERATION %d/%d] <<<\n", iter, max_iterations))
  cat(sprintf("    Testing theta: [%.4f, %.4f, %.4f, %.4f, %.4f]\n",
              current_theta[1], current_theta[2], current_theta[3],
              current_theta[4], current_theta[5]))
  
  pct_below_threshold_e_sum = 0
  pct_below_threshold_p_sum = 0
  
  for (reps in 1:num_reps){
    # Reset simulation
    reset_simulation_state()
    
    source('./MISC/MAIN_SIM_LOOP_OVER_T_DIFFUSION.R')
    
    longitudinal_df = data.frame(
      epithelium_longitudinal,
      macrophages_longitudinal,
      microbes_longitudinal,
      tregs_longitudinal,
      microbes_cumdeath_longitudinal,
      pathogens_lumen_longitudinal
    )
    
    colnames(longitudinal_df) = c(colnames_insert,'pathogens_lumen')
    
    longitudinal_df$t = 1:t_max
    longitudinal_df$sterile = sterile
    longitudinal_df$macspec_on = macspec_on
    longitudinal_df$tregs_on = allow_tregs
    longitudinal_df$ros_level = ros_level
    longitudinal_df$pat_level = pat_level
    longitudinal_df$randomize_tregs = randomize_tregs
    longitudinal_df$param_set_id = param_set_use$param_set_id
    longitudinal_df$rep_id = reps
    
    current_epithelial_score = longitudinal_df$epithelial_score
    current_pathogens        = longitudinal_df$pathogen
    
    sum_e = sum(round(current_epithelial_score, 3))
    sum_p = sum(round(current_pathogens, 3))
    
    recent_scores_success_e = tail(current_epithelial_score, success_duration)
    recent_scores_success_p = tail(current_pathogens, success_duration)
    
    recent_scores_success_e = round(recent_scores_success_e, 3)
    recent_scores_success_p = round(recent_scores_success_p, 3)
    
    pct_above_threshold_e = round(sum(recent_scores_success_e < success_threshold_e) / length(recent_scores_success_e),3)
    pct_below_threshold_p = round(sum(recent_scores_success_p < success_threshold_p) / length(recent_scores_success_p),3)
    
    is_success_e = (pct_above_threshold_e > success_rate)
    is_success_p = (pct_below_threshold_p > success_rate)
    
    cat(paste0(param_set_id_use," ", reps, " ", pat_level, " ", ros_level, " ", overwrite_in, " "),file = success_log_file, append = TRUE)
    cat(paste0(current_theta[1]," ", 
               current_theta[2]," ",
               current_theta[3]," ", 
               current_theta[4]," ",
               current_theta[5]," "),file = success_log_file, append = TRUE)
    cat(paste0(round(min(recent_scores_success_e),3)," ", 
               round(min(recent_scores_success_p),3)," ",
               round(mean(recent_scores_success_e),3)," ", 
               round(mean(recent_scores_success_p),3)," ",
               pct_above_threshold_e," ",
               pct_below_threshold_p," ",
               sum_e," ",
               sum_p," ",
               "\n"),file = success_log_file, append = TRUE)
    
    pct_below_threshold_e_sum = pct_below_threshold_e_sum + pct_above_threshold_e
    pct_below_threshold_p_sum = pct_below_threshold_p_sum + pct_below_threshold_p
    
    # if(reps==2 & pct_below_threshold_e_sum==0){
    #   print('this one not gonna happen')
    #   break; #stop trying
    # }

    
    # saveRDS(current_theta, paste0(dir_name_data,
    #                               '/theta_param_set_id_',param_set_id_use,
    #                               '_overwrite_',overwrite_in,
    #                               '_sterile_',sterile,
    #                               '_macspec_',macspec_on,
    #                               '_tregs_',allow_tregs,
    #                               '_ros_level_',ros_level,
    #                               '_pat_level_',pat_level,
    #                               '_trnd_',randomize_tregs,
    #                               '_n1_',n1,
    #                               '_n2_',n2,
    #                               '_iter_',iter,
    #                               '_rep_',reps,
    #                               '.rds'))
    
    # ## SAVE TO COME BACK TO IT LATER?
    # saveRDS(longitudinal_df, paste0(dir_name_data,
    #                                 '/longitudinal_df_param_set_id_',param_set_id_use,
    #                                 '_overwrite_',overwrite_in,
    #                                 '_sterile_',sterile,
    #                                 '_macspec_',macspec_on,
    #                                 '_tregs_',allow_tregs,
    #                                 '_ros_level_',ros_level,
    #                                 '_pat_level_',pat_level,
    #                                 '_trnd_',randomize_tregs,
    #                                 '_n1_',n1,
    #                                 '_n2_',n2,
    #                                 '_iter_',iter,
    #                                 '_rep_',reps,
    #                                 '.rds'))
    
  }
}
