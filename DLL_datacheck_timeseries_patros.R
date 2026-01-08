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

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

save_images_path = 'timeseries_tri_all_param_ids'
dir.create(paste0("/Users/burcutepekule/Desktop/",save_images_path), showWarnings = FALSE)
# path  = "/Users/burcutepekule/Desktop/sim_abm/"
path  = "/Users/burcutepekule/Desktop/sim_abm_local/"

# Define the ros and pat value ranges
ros_vals        = c(1,3,5,10) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
# ros_vals        = c(0,1,3,5,10) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
pat_vals        = c(1,2,5,10,100)
tregs_on_in     = 1
sterile_in      = 0

# File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
# Dynamically create file lists for all ros/pat combinations
all_files = list()
all_indices = list()
for (ros in ros_vals) {
  for (pat in pat_vals) {
    var_name = paste0("files_", ros, "_", pat)
    pattern_str = paste0("^longitudinal_df_param_set_id_\\d+\\_sterile_0_macspec_0_tregs_",tregs_on_in,"_ros_level_", ros, "_pat_level_", pat, "_trnd_0.rds$")
    files_temp = list.files(path, pattern = pattern_str, full.names = TRUE)
    all_files[[var_name]] = files_temp
    
    # Extract indices
    indices_name = paste0("indices_", ros, "_", pat)
    all_indices[[indices_name]] = str_extract(basename(files_temp), "\\d+") |> as.numeric()
  }
}

# Find common indices across all combinations
# indices = Reduce(intersect, all_indices)
indices = 50002

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

# ===== reread to include local results 
param_id_vec = indices
rep_ind_vec  = 0:2 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
alpha_plot   = 1/length(rep_ind_vec)
i_opt        = 1
# variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
# variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2"),c("phagocyte_M1M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
variables_2_plot = list("epithelial_score","pathogen")
background_on    = c(1,1,rep(0,length(variables_2_plot)-2))

epithelial_limit = 5
pathogen_limit   = 10
max_level_injury = 10*25 # 2*x0_in*grid_size
triangle_df = c()

for(param_id in param_id_vec){
  
  source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
  TF_matrix = control_matrix
  
  # last_row = control_matrix[dim(control_matrix)[1],]
  # triangle_df = rbind(triangle_df, c(param_id, is_inverted_triangle_boundary(control_matrix), is_rectangular_top(control_matrix), all(last_row==FALSE)))
  
  TF_df = as.data.frame.table(TF_matrix, stringsAsFactors = FALSE)
  colnames(TF_df) = c("pat_level", "ros_level", "controlled")
  
  # Extract numeric values from row/column names
  TF_df$pat_level = as.numeric(gsub("pat", "", TF_df$pat_level))
  TF_df$ros_level = as.numeric(gsub("ros", "", TF_df$ros_level))
  
  # Combine all results
  results = do.call(rbind, results_list)
  
  results = results %>% dplyr::filter(rep_id %in% rep_ind_vec)
  results = results %>% dplyr::mutate(phagocyte_M1M2 = phagocyte_M1+phagocyte_M2)
  
  for (p_ind in 1:length(variables_2_plot)){
    variables = variables_2_plot[p_ind][[1]]
    
    data_long = results %>%
      dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    if(background_on[p_ind]==1){
      p = ggplot(data_long, aes(x = t, y = value)) +
        # Background
        geom_rect(data = TF_df,
                  aes(fill = controlled),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                  alpha = 0.1, inherit.aes = FALSE) +
        scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red"),
                          name = "Under Control") +
        # Border
        geom_rect(data = TF_df,
                  aes(color = controlled),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                  fill = NA, linewidth = 1.5, inherit.aes = FALSE,
                  show.legend = FALSE) +
        scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkred")) +
        # Reset color scale for lines
        new_scale_color() +
        # Horizontal threshold line
        geom_hline(yintercept = max_level_injury,
                   linetype = "solid",
                   color = "black",
                   linewidth = 0.5) +
        geom_hline(yintercept = epithelial_limit,
                   linetype = "solid",
                   color = "black",
                   linewidth = 0.5) +
        # Lines with agent colors
        geom_line(aes(color = variable, group = rep_id), 
                  alpha = alpha_plot, linewidth = 1) +
        scale_color_manual(values = agent_colors, name = "Agent") +
        facet_grid(pat_level ~ ros_level, labeller = label_both) +
        theme_minimal() +
        scale_y_log10(limits = c(NA, max_level_injury)) +
        # labs(title = title_opt, x = "Time", y = "Count")
        labs(title = "", x = "Time", y = "Count")
    }else{
      p = ggplot(data_long, aes(x = t, y = value))+
        geom_hline(yintercept = 10,
                   linetype = "solid",
                   color = "gray",
                   linewidth = 0.5) +
        # Lines with agent colors
        geom_line(aes(color = variable, group = interaction(rep_id, variable)), 
                  alpha = alpha_plot, linewidth = 1) +
        scale_color_manual(values = agent_colors, name = "Agent") +
        facet_grid(pat_level ~ ros_level, labeller = label_both) +
        theme_minimal() +
        # labs(title = title_opt, x = "Time", y = "Count")
        labs(title = "", x = "Time", y = "Count")
    }
    
    ggsave(
      # filename = paste0("/Users/burcutepekule/Desktop/",save_images_path,"/sterile_",sterile_in,"_tregs_on_",tregs_on_in,"_",param_id,"_",variables[1],".png"),
      filename = paste0("/Users/burcutepekule/Desktop/",save_images_path,"/i_opt_",i_opt,"_sterile_",sterile_in,"_tregs_on_",tregs_on_in,"_",param_id,"_",variables[1],".png"),
      plot = p,
      width = 24,
      height = 16,
      dpi = 300,
      bg='white'
    )
  }
} 


# saveRDS(triangle_df, 'triangle_df.rds')
# 
# triangle_df = as.data.frame(readRDS('triangle_df.rds'))
# colnames(triangle_df) = c('param_id','is_triangle','is_last_row_all_false','is_rectangular')
# 
# triangle_df_1 = triangle_df %>% dplyr::filter((is_triangle & is_last_row_all_false)|(is_rectangular & is_last_row_all_false))
