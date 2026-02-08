rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggnewscale)
library(stringr)
library(zoo)
library(philentropy)
library(purrr)
library(data.table)

setwd('~/Dropbox/tregs_clean/')

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")
source("./MISC/LOAD_PAT_LEVELS_DFF.R") # loads pat_level_vectors

indices_vec = c(147, 172, 222, 269)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
split_equal = function(x, n_chunks) {
  if (n_chunks == 1) {
    return(list(`1` = x))
  }
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

save_data_path_data = "/Users/burcutepekule/Desktop/ts_diffusion_highres_eff_data"

dir.create(paste0(save_data_path_data)) # for desktop

epithelial_limit      = 5
pathogen_limit        = 0.1 # max is 1
max_level_injury      = 10*25 # 2*x0_in*grid_size - THIS STILL HOLDS FOR DIFFUSION!

trnd_in      = 0
sterile_in   = 0
macspec_in   = 0
overwrite_in = 0
rep_ind_vec  = 0:9 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
pat_vals     = seq(1,5,0.5)

scenarios_df = c()
### Tregs OFF case
scenarios_df = rbind(scenarios_df, expand.grid(
  param_id    = indices_vec,
  tregs_on_in = 0,
  i_opt       = 0,
  treg_eff_in = 0,
  m2on_in     = 0
))
### Tregs ON case
scenarios_df = rbind(scenarios_df, expand.grid(
  param_id    = indices_vec,
  tregs_on_in = 1,
  i_opt       = 1:5,
  # treg_eff_in = 0:10,
  # m2on_in     = c(0,1)
  treg_eff_in = 10,
  m2on_in     = 1
))
dim(scenarios_df)[1]

chunks       = split_equal(1:nrow(scenarios_df), n1)
loop_over_sc = chunks[[n2]]

for (scenario_ind in loop_over_sc){
  
  param_id    = scenarios_df[scenario_ind,]$param_id
  tregs_on_in = scenarios_df[scenario_ind,]$tregs_on_in
  i_opt       = scenarios_df[scenario_ind,]$i_opt
  treg_eff_in = scenarios_df[scenario_ind,]$treg_eff_in
  m2on_in     = scenarios_df[scenario_ind,]$m2on_in
  
  path_data  = paste0("/Users/burcutepekule/Desktop/sim_dff_highres_",param_id,"/")
  
  col_avg_keep = c()
  ros_vals     = ros_level_vectors[[as.character(param_id)]] # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  
  message(paste0("Processing param_id ", param_id, " for: ", 
                 " tregs_on_in: ", tregs_on_in, 
                 " i_opt: ", i_opt, 
                 " treg_eff_in: ", treg_eff_in, 
                 " m2on_in: ", m2on_in))
  
  source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros_diff_eff.R')
  
  control_matrix_long = as.data.frame(control_matrix_long)
  colnames(control_matrix_long) = c('param_id','pat','ros','pct')
  TF_matrix = control_matrix
  
  TF_df = as.data.frame.table(TF_matrix, stringsAsFactors = FALSE)
  colnames(TF_df) = c("pat_level", "ros_level", "control_percentage")
  
  # Extract numeric values from row/column names
  TF_df$pat_level = as.numeric(gsub("pat", "", TF_df$pat_level))
  TF_df$ros_level = as.numeric(gsub("ros", "", TF_df$ros_level))
  
  # Convert control_percentage to numeric (it's a proportion from 0 to 1)
  TF_df$control_percentage = as.numeric(TF_df$control_percentage)
  
  # Create a label showing the percentage
  TF_df$control_label = paste0(round(TF_df$control_percentage * 100), "%")
  
  # Calculate average percentage for each column (ros_level)
  col_avg = TF_df %>%
    group_by(ros_level) %>%
    # summarise(avg_pct = mean(control_percentage, na.rm = TRUE)) %>%
    summarise(avg_pct = sum(pat_level*control_percentage, na.rm = TRUE)/sum(pat_vals)) 
  
  col_avg$opt_ind     = i_opt
  col_avg$overwrite   = overwrite_in
  col_avg$treg_eff    = treg_eff_in
  
  results = results %>% dplyr::filter(rep_id %in% rep_ind_vec)
  
  save(TF_matrix, col_avg, results, file = paste0(save_data_path_data, 
                                                  "/img_",param_id,
                                                  "_i_opt_",i_opt,
                                                  "_effidx_",treg_eff_in,
                                                  "_m2on_",m2on_in,
                                                  ".RData"))
  
  message(paste0("SAVED param_id ", param_id, " for: ", 
                 " tregs_on_in: ", tregs_on_in, 
                 " i_opt: ", i_opt, 
                 " treg_eff_in: ", treg_eff_in, 
                 " m2on_in: ", m2on_in))
  
}

