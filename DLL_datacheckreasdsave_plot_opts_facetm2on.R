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
path_data  = "/Users/burcutepekule/Desktop/data_diffusion"

variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2"),"treg_active")
background_on    = c(1,1,rep(0,length(variables_2_plot)-2))
epithelial_limit = 5
pathogen_limit   = 0.1 # max is 1
max_level_injury = 10*25 # 2*x0_in*grid_size - THIS STILL HOLDS FOR DIFFUSION!
alpha_plot       = 0.2

for (param_id in c(147, 172, 222, 269)){
  
  print(param_id)
  path_img = paste0("/Users/burcutepekule/Desktop/ts_diffusion_highres_",param_id)
  dir.create(path_img, showWarnings = TRUE)
  
  # first load 0 - add to both facets
  col_avg_keep = c()
  load(paste0(path_data,'/img_',param_id,'_i_opt_',0,'_effidx_',0,'_m2on_',0,'.RData'))
  col_avg$m2on = 0
  col_avg_keep = rbind(col_avg_keep, col_avg)
  col_avg$m2on = 1
  col_avg_keep = rbind(col_avg_keep, col_avg)
  
  # then opts
  effidx    = 10
  i_opt_vec = 1
  m2on_vec  = c(0,1)
  for (m2on in m2on_vec){
    for(i_opt in i_opt_vec){
      print(i_opt)
      load(paste0(path_data,'/img_',param_id,'_i_opt_',i_opt,'_effidx_',effidx,'_m2on_',m2on,'.RData'))
      col_avg$m2on = m2on
      col_avg_keep = rbind(col_avg_keep, col_avg)
      gc()
    }
  }
  
  # Create color palette: red for 0, viridis for the rest
  n_other    = length(i_opt_vec)   # number of non-zero values
  blues_func = colorRampPalette(RColorBrewer::brewer.pal(n_other + 1, "Blues"))
  blues_cols = blues_func(n_other + 1)[-1]
  custom_colors = c("0" = "red", setNames(blues_cols, as.character(i_opt_vec)))
  
  p_opt = ggplot(col_avg_keep, aes(x = ros_level, y = avg_pct, color = factor(opt_ind))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(
      trans = "exp",  # This creates exponential spacing
      breaks = unique(col_avg_keep$ros_level)  # Use actual ros_level values as breaks
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "",
      x = "ROS Level",
      y = "Average % Controlled",
      color = "opt_ind"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      strip.text = element_text(size = 12)
    ) + 
    scale_color_manual(values = custom_colors)+
    facet_wrap(~ m2on, labeller = label_both, nrow = 1)
  
  ggsave(
    filename = paste0(path_img, "/param_",param_id,"_OPTS_m2on.png"),
    plot = p_opt,
    width = 16,
    height = 8,
    dpi = 300,
    bg='white'
  )
  
}
