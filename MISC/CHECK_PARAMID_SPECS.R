rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# CHUNKS
# ============================================================================

# ============================================================================
loop_over = c(81250, 88750, 30000, 45500, 92000, 50250, 62500, 68752)
# ============================================================================

params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)

params_df[order(params_df$diffusion_speed_DAMPs),c('param_set_id','diffusion_speed_DAMPs')]
params_df[order(params_df$diffusion_speed_ROS),c('param_set_id','diffusion_speed_ROS')]
params_df[order(params_df$DAMPs_decay),c('param_set_id','DAMPs_decay')]
params_df[order(params_df$SAMPs_decay),c('param_set_id','SAMPs_decay')]
params_df[order(params_df$PAMPs_decay),c('param_set_id','PAMPs_decay')]
params_df[order(params_df$recruitment_rate_danger),c('param_set_id','recruitment_rate_danger')]


params_df[order(params_df$PAMPs_decay),c('param_set_id','PAMPs_decay','recruitment_rate_danger',
                                         'diffusion_speed_DAMPs','SAMPs_decay')]


params_df_long = pivot_longer(params_df, cols = -c('param_set_id'), names_to='parameter', values_to = 'value')
