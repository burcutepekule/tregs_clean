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
path_data  = "/Users/burcutepekule/Desktop/data_diffusion_rpat05*leak/"

trnd_in      = 0
sterile_in   = 0
macspec_in   = 0
overwrite_in = 0
rep_ind_vec  = 0:9 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
pat_vals     = seq(1,5,0.5)

variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2"),"treg_active")
background_on    = c(1,1,rep(0,length(variables_2_plot)-2))
source('./MISC/PERFORMANCE_METRICS.R')
alpha_plot = 1

param_id    = 147
# tregs_on_in = 1
# i_opt       = 5
# treg_eff_in = 10
# m2on_in     = 1

### tregs off
tregs_on_in = 0
i_opt       = 0
treg_eff_in = 0
m2on_in     = 0

# pat_level_1box = 2
# ros_level_1box = 3.9

pat_level_1box = 5
ros_level_1box = 7.8

message(paste0("Printing grid for param_id ", param_id, " for: ", 
               " tregs_on_in: ", tregs_on_in, 
               " i_opt: ", i_opt, 
               " treg_eff_in: ", treg_eff_in, 
               " m2on_in: ", m2on_in))


path_img = paste0("/Users/burcutepekule/Desktop/ts_diffusion_highres_",param_id)
dir.create(path_img, showWarnings = TRUE)

load(paste0(path_data,'/img_',param_id,'_i_opt_',i_opt,'_effidx_',treg_eff_in,'_m2on_',m2on_in,'.RData'))

results$ros_level = round(results$ros_level,3)
results = results %>% dplyr::filter(ros_level==ros_level_1box & pat_level==pat_level_1box)

for (p_ind in 1:length(variables_2_plot)){
  source('./MISC/SUBPLOT_ONE_GRID.R')
  ggsave(
    filename = paste0(path_img,
                      "/param_",param_id,
                      "_i_opt_",i_opt,
                      '_effidx_',treg_eff_in,
                      '_m2on_',m2on_in,
                      "_",variables_name,"_1box.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    bg='white'
  )
}


results_r = results %>% dplyr::filter(rep_id==1)

