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

# 45500, 92000, 30000, 81250, 88750, 50250
indices = c(45500, 92000, 30000, 81250, 88750, 50250, 62500)## you can either do this at the very top or pick the indices that are complete with the loop below

save_images_path_data = "timeseries_tri_all_param_ids_suppress"
dir.create(paste0("/Users/burcutepekule/Desktop/",save_images_path_data))

# path_data  = "/Users/burcutepekule/Desktop/sim_abm_50250_88750/"
path_data  = "/Users/burcutepekule/Desktop/sim_abm/"

# Define the ros and pat value ranges
ros_vals        = seq(0,10,1) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
trnd_in         = 0
sterile_in      = 0
macspec_in      = 0
overwrite_in_vec= c(0,1)

# source("~/Dropbox/tregs_clean/MISC/FIND_COMPLETE_PARAMIDS.R")

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

# ===== reread to include local results
rep_ind_vec  = 0:9 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
alpha_plot   = 2/length(rep_ind_vec)

# Function to dynamically determine i_opt_vec based on param_id
get_i_opt_vec = function(param_id, path_data) {
  # Default mapping based on known param_id optimization ranges
  i_opt_map = list(
    "50250" = 0:4,   # up to 4 for 50250
    "88750" = 0:8,   # up to 8 for 88750
    "62500" = 0:4,   # similar parameter set, assume similar range
    "default" = 0    # 0 means no optimization for others
  )

  # Check if param_id has a specific mapping
  param_id_str = as.character(param_id)
  if (param_id_str %in% names(i_opt_map)) {
    return(i_opt_map[[param_id_str]])
  }

  # Otherwise, try to detect from available files in path_data
  if (dir.exists(path_data)) {
    files = list.files(path_data, pattern = paste0("longitudinal_df_param_set_id_", param_id, "_.*_optidx_"))
    if (length(files) > 0) {
      # Extract optidx values from filenames
      opt_indices = unique(as.numeric(gsub(".*_optidx_(\\d+)\\.rds$", "\\1", files)))
      opt_indices = opt_indices[!is.na(opt_indices)]
      if (length(opt_indices) > 0) {
        return(sort(opt_indices))
      }
    }
  }

  # Fallback to default
  return(i_opt_map[["default"]])
}

# variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
# background_on    = c(1,1,rep(0,length(variables_2_plot)-2))

variables_2_plot      = list("epithelial_score")
background_on         = c(1)
epithelial_limit      = 5
pathogen_limit        = 10
max_level_injury      = 10*25 # 2*x0_in*grid_size
triangle_df           = c()

# if(file.exists('./ids_read_abm_patros.rds')){
#   inds_read = readRDS('./ids_read_abm_patros.rds')
# }else{
#   inds_read = c()
# }
inds_read = c()
inds2read = sort(setdiff(indices,inds_read))
length(inds2read)
processed_indices = c()
plot_in = 1

for(param_id in inds2read){

  if(param_id %in% c(30000, 81250)){
    pat_vals = c(1, 2, seq(2.5, 5, 0.5))
  }else if(param_id %in% c(45500, 88750, 50250, 62500)){
    pat_vals = c(1, 2, seq(5, 10, 1))
  }else if(param_id %in% c(92000)){
    pat_vals = c(1, 5, 7,seq(8, 12, 1))
  }

  # Dynamically determine i_opt_vec based on param_id
  i_opt_vec = get_i_opt_vec(param_id, path_data)
  message(paste0("Processing param_id ", param_id, " with i_opt_vec: ", paste(i_opt_vec, collapse = ", ")))

  for(i_opt in i_opt_vec){
    
    if(i_opt ==0){
      tregs_on_in = 0
    }else{
      tregs_on_in = 1
    }
    
    for(overwrite_in in overwrite_in_vec){
      source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
      
      control_matrix_long = as.data.frame(control_matrix_long)
      colnames(control_matrix_long) = c('param_id','pat','ros','pct')
      
      processed_indices = c(processed_indices, param_id) #add 
      TF_matrix = control_matrix
      
      control_matrix_bin = control_matrix>0.5
      last_row = control_matrix_bin[dim(control_matrix_bin)[1],]
      
      triangle_df = rbind(triangle_df, c(param_id, is_inverted_triangle_boundary(control_matrix_bin), is_rectangular_top(control_matrix_bin), all(last_row==FALSE)))
      
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
      
      if(plot_in==1){
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
            filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,
                              "/param_",param_id,
                              "_i_opt_",i_opt,
                              "_overwrite_",overwrite_in,
                              "_sterile_", sterile_in,
                              "_tregs_on_",tregs_on_in,
                              "_tregs_rnd_",trnd_in,
                              "_",variables[1],".png"),
            plot = p,
            width = 24,
            height = 16,
            dpi = 300,
            bg='white'
          )
        }}
    } 
  }
}


# saveRDS(triangle_df, 'triangle_df.rds')
# # #   
# triangle_df = as.data.frame(readRDS('triangle_df.rds'))
# colnames(triangle_df) = c('param_id','is_triangle','is_rectangular','is_last_row_all_false')
# 
# triangle_df_1 = triangle_df %>% dplyr::filter((is_triangle)|(is_rectangular & is_last_row_all_false))
# 
# ### pick the one with the smallest active age limit?
# 
# params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
# params_df_1 = params_df %>% dplyr::filter(param_set_id %in% triangle_df_1$param_id)
# params_df_1 = params_df_1[c('param_set_id','active_age_limit')]
# params_df_1 = params_df_1 %>% arrange(active_age_limit)
# params_df_1 = params_df_1 %>% arrange(param_set_id)
# 
# # Update the list of read indices
# if(!file.exists('./ids_read_abm_patros.rds')){
#   saveRDS(processed_indices, './ids_read_abm_patros.rds')
# }else{
#   inds_read_old = readRDS('./ids_read_abm_patros.rds')
#   inds_read_updated = c(inds_read_old, processed_indices)
#   saveRDS(inds_read_updated, './ids_read_abm_patros.rds')
# }
