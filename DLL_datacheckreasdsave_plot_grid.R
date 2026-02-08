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
path_data  = "/Users/burcutepekule/Desktop/ts_diffusion_highres_eff_data/"

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

variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2"),"treg_active")
background_on    = c(1,1,rep(0,length(variables_2_plot)-2))
epithelial_limit = 5
pathogen_limit   = 0.1 # max is 1
max_level_injury = 10*25 # 2*x0_in*grid_size - THIS STILL HOLDS FOR DIFFUSION!
alpha_plot       = 0.2

for (scenario_ind in loop_over_sc){

  param_id    = scenarios_df[scenario_ind,]$param_id
  tregs_on_in = scenarios_df[scenario_ind,]$tregs_on_in
  i_opt       = scenarios_df[scenario_ind,]$i_opt
  treg_eff_in = scenarios_df[scenario_ind,]$treg_eff_in
  m2on_in     = scenarios_df[scenario_ind,]$m2on_in
  
  
  message(paste0("Printing grid for param_id ", param_id, " for: ", 
                 " tregs_on_in: ", tregs_on_in, 
                 " i_opt: ", i_opt, 
                 " treg_eff_in: ", treg_eff_in, 
                 " m2on_in: ", m2on_in))
  
  
  path_img = paste0("/Users/burcutepekule/Desktop/ts_diffusion_highres_",param_id)
  dir.create(path_img, showWarnings = TRUE)
  
  load(paste0(path_data,'/img_',param_id,'_i_opt_',i_opt,'_effidx_',treg_eff_in,'_m2on_',m2on_in,'.RData'))
  
  for (p_ind in 1:length(variables_2_plot)){
    source('./MISC/SUBPLOT_ONE_GRID.R')
    ggsave(
      filename = paste0(path_img,
                        "/param_",param_id,
                        "_i_opt_",i_opt,
                        '_effidx_',treg_eff_in,
                        '_m2on_',m2on_in,
                        "_",variables_name,".png"),
      plot = p,
      width = 48,
      height = 24,
      dpi = 300,
      bg='white'
    )
  }
  
}
