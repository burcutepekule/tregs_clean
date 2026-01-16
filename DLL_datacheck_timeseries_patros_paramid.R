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

setwd('~/Dropbox/tregs_clean/')

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

indices          = 80750 # 62500, 73750, 80750 - 67250
save_images_path_data = paste0('timeseries_tri_all_param_ids_',indices)
dir.create(paste0("/Users/burcutepekule/Desktop/",save_images_path_data), showWarnings = FALSE)

path_data  = "/Users/burcutepekule/Desktop/sim_abm/"

# Define the ros and pat value ranges
ros_vals        = seq(0,10,1) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
pat_vals        = c(1, 2, 5, 6, 7, 8, 9, 10)
# pat_vals        = c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5)
tregs_on_in     = 1
trnd_in         = 0
sterile_in      = 0
macspec_in      = 0

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

# ===== reread to include local results 
rep_ind_vec  = 0:4 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
alpha_plot   = 2/length(rep_ind_vec)
i_opt        = NA

variables_2_plot = list("epithelial_score")
background_on    = c(1)

epithelial_limit      = 5
pathogen_limit        = 10
max_level_injury      = 10*25 # 2*x0_in*grid_size
triangle_df           = c()

for(param_id in indices){
  
  source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
  TF_matrix = control_matrix

  last_row = control_matrix[dim(control_matrix)[1],]
  # triangle_df = rbind(triangle_df, c(param_id, is_inverted_triangle_boundary(control_matrix), is_rectangular_top(control_matrix), all(last_row==FALSE)))

  TF_df = as.data.frame.table(TF_matrix, stringsAsFactors = FALSE)
  colnames(TF_df) = c("pat_level", "ros_level", "control_percentage")

  # Extract numeric values from row/column names
  TF_df$pat_level = as.numeric(gsub("pat", "", TF_df$pat_level))
  TF_df$ros_level = as.numeric(gsub("ros", "", TF_df$ros_level))

  # Convert control_percentage to numeric (it's a proportion from 0 to 1)
  TF_df$control_percentage = as.numeric(TF_df$control_percentage)

  # Create a label showing the percentage
  TF_df$control_label = paste0(round(TF_df$control_percentage * 100), "%")

  # Combine all results
  results = do.call(rbind, results_list)
  
  results = results %>% dplyr::filter(rep_id %in% rep_ind_vec)
  # results = results %>% dplyr::mutate(phagocyte_M1M2 = phagocyte_M1+phagocyte_M2)
  
  for (p_ind in 1:length(variables_2_plot)){
    variables = variables_2_plot[p_ind][[1]]
    
    data_long = results %>%
      dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    if(background_on[p_ind]==1){
      p = ggplot(data_long, aes(x = t, y = value)) +
        # Background with gradient coloring
        geom_rect(data = TF_df,
                  aes(fill = control_percentage),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                  alpha = 0.3, inherit.aes = FALSE) +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                            midpoint = 0.5,
                            limits = c(0, 1),
                            name = "% Controlled",
                            labels = scales::percent) +
        # Border with gradient coloring
        geom_rect(data = TF_df,
                  aes(color = control_percentage),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                  fill = NA, linewidth = 1.5, inherit.aes = FALSE,
                  show.legend = FALSE) +
        scale_color_gradient2(low = "darkred", mid = "orange", high = "darkgreen",
                             midpoint = 0.5,
                             limits = c(0, 1)) +
        # Add percentage text labels in top right corner of each facet
        geom_text(data = TF_df,
                  aes(label = control_label),
                  x = Inf, y = Inf,
                  hjust = 1.1, vjust = 1.5,
                  size = 4, fontface = "bold",
                  color = "black",
                  inherit.aes = FALSE) +
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
        scale_y_log10() +
        # scale_y_log10(limits = c(NA, max_level_injury)) +
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
      # filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,"/sterile_",sterile_in,"_tregs_on_",tregs_on_in,"_",param_id,"_",variables[1],".png"),
      filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,"/i_opt_",i_opt,
                        "_sterile_", sterile_in,
                        "_tregs_on_",tregs_on_in,
                        "_tregs_rnd_",trnd_in,
                        "_",param_id,"_",variables[1],".png"),
      plot = p,
      width = 24,
      height = 16,
      dpi = 300,
      bg='white'
    )
  }
} 