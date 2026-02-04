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
source("./MISC/LOAD_PAT_LEVELS.R") # loads pat_level_vectors
# path_data  = "/Users/burcutepekule/Desktop/sim_abm/"
path_data  = "/Users/burcutepekule/Desktop/sim_diffusion_local/"

# indices_vec = c(43596)
# indices_vec = c(43597:43606)
indices_vec = c(43601:43606)
# indices_vec = c(43591, 43604)
indices_vec = c(43604)
indices_vec = c(43606:43676) 
indices_vec = c(43672)
# indices_vec = 43600 + c(30, 46)
indices_vec = 43600 + c(46)


# ros_vals = seq(0,10,1)
# pat_vals = c(1,5,10,15,20)
# # File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
# # Dynamically create file lists for all ros/pat combinations
# all_files = list()
# all_indices = list()
# for (ros_in in ros_vals) {
#   for (pat_in in pat_vals) {
#     var_name = paste0("files_", ros_in, "_", pat_in)
#     pattern_str = paste0("^longitudinal_df_param_set_id_\\d+\\_sterile_",0,
#                          "_macspec_",0,
#                          "_tregs_",0,
#                          "_ros_level_",ros_in,
#                          "_pat_level_",pat_in,
#                          "_trnd_",0,
#                          "_overwrite_",0,
#                          "_optidx_",0,
#                          ".rds$")
#     files_temp = list.files(path_data, pattern = pattern_str, full.names = TRUE)
#     all_files[[var_name]] = files_temp
#     
#     # Extract indices
#     indices_name = paste0("indices_", ros_in, "_", pat_in)
#     all_indices[[indices_name]] = str_extract(basename(files_temp), "\\d+") |> as.numeric()
#   }
# }
# 
# # Find common indices across all combinations
# # indices = Reduce(intersect, all_indices)
# indices_vec = sort(Reduce(union, all_indices))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
split_equal = function(x, n_chunks) {
  if (n_chunks == 1) {
    return(list(`1` = x))
  }
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

# args   = commandArgs(trailingOnly = TRUE)
# n1     = as.integer(args[1])
# n2     = as.integer(args[2])
n1     = 1
n2     = 1

chunks       = split_equal(indices_vec, n1)
loop_over_sc = chunks[[n2]]

skip_opt_0  = 0
save_images_path_data = "ts_diffusion"
dir.create(paste0("/Users/burcutepekule/Desktop/",save_images_path_data))

# path_data  = "/Users/burcutepekule/Desktop/sim_abm_local/"

for (indices in loop_over_sc){
  # Define the ros and pat value ranges
  ros_vals        = seq(0,10,1) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
  trnd_in         = 0
  sterile_in      = 0
  macspec_in      = 0
  overwrite_in_vec= c(0, 1)
  
  # source("~/Dropbox/tregs_clean/MISC/FIND_COMPLETE_PARAMIDS.R")
  # ===== reread to include local results
  rep_ind_vec  = 0:9 # this limits the max rep index read by the data in source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
  alpha_plot   = 2/length(rep_ind_vec)
  
  # variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
  
  variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2","phagocyte_M0"))
  background_on    = c(1,1,rep(0,length(variables_2_plot)-2))
  
  # variables_2_plot      = list("epithelial_score")
  # background_on         = c(1)
  
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
  skip_id = 0
  
  for(param_id in inds2read){
    
    col_avg_keep = c()
    pat_level_vec   = pat_level_vectors[[as.character(param_id)]]
    
    # Dynamically determine pat_vals and i_opt_vec based on param_id
    pat_vals = get_pat_vals(param_id, path_data) # from /MISC/PLOT_FUNCTIONS_ABM.R
    i_opt_vec = get_i_opt_vec(param_id, path_data) # from /MISC/PLOT_FUNCTIONS_ABM.R

    if(length(i_opt_vec)<1){
      message("Skipping id ", param_id, " due to data")
      skip_id = 1
    }else{
      
      if(skip_opt_0==1){
        i_opt_vec = i_opt_vec[i_opt_vec>0] # IF YOU ALREADY PRINTED THE VANILLA!
      }
      
      message(paste0("Processing param_id ", param_id, " with pat_vals: ", paste(pat_vals, collapse = ", "), 
                     " and i_opt_vec: ", paste(i_opt_vec, collapse = ", ")))
      
      for(i_opt in i_opt_vec){
        
        if(i_opt ==0){
          tregs_on_in = 0
          prefix      = '0_'
        }else{
          tregs_on_in = 1
          prefix      = ''
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
          
          # Calculate average percentage for each column (ros_level)
          col_avg = TF_df %>%
            group_by(ros_level) %>%
            # summarise(avg_pct = mean(control_percentage, na.rm = TRUE)) %>%
            summarise(avg_pct = sum(pat_level*control_percentage, na.rm = TRUE)/sum(pat_level_vec)) 
          
          col_avg$opt_ind     = i_opt
          col_avg$overwrite   = overwrite_in
          
          col_avg_keep = rbind(col_avg_keep, col_avg)
          
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
              
              # Create custom labeller that includes column averages
              col_labeller = function(labels) {
                lapply(labels, function(x) {
                  avg = col_avg$avg_pct[col_avg$ros_level == x]
                  paste0("ros_level: ", x, "\n(avg: ", scales::percent(avg, accuracy = 0.1), ")")
                })
              }
              
              if(background_on[p_ind]==1){
                if(variables=='epithelial_score'){
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
                    # facet_grid(pat_level ~ ros_level, labeller = label_both) +
                    # replace labeller to show column averages
                    facet_grid(pat_level ~ ros_level, 
                               labeller = labeller(pat_level = label_both, 
                                                   ros_level = col_labeller))+
                    theme_minimal() +
                    scale_y_log10() +
                    # scale_y_log10(limits = c(NA, max_level_injury)) +
                    # labs(title = title_opt, x = "Time", y = "Count")
                    labs(title = "", x = "Time", y = "Count")+
                    # Increase axis text size
                    theme(
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
                      strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
                      legend.text = element_text(size = 14),   # legend item labels
                      legend.title = element_text(size = 16)   # legend title
                    )
                }else if(variables=='pathogen'){
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
                    geom_hline(yintercept = pathogen_limit,
                               linetype = "solid",
                               color = "black",
                               linewidth = 0.5) +
                    # Lines with agent colors
                    geom_line(aes(color = variable, group = rep_id),
                              alpha = alpha_plot, linewidth = 1) +
                    scale_color_manual(values = agent_colors, name = "Agent") +
                    # facet_grid(pat_level ~ ros_level, labeller = label_both) +
                    # replace labeller to show column averages
                    facet_grid(pat_level ~ ros_level, 
                               labeller = labeller(pat_level = label_both, 
                                                   ros_level = col_labeller))+
                    theme_minimal() +
                    scale_y_log10() +
                    # scale_y_log10(limits = c(NA, max_level_injury)) +
                    # labs(title = title_opt, x = "Time", y = "Count")
                    labs(title = "", x = "Time", y = "Count")+
                    # Increase axis text size
                    theme(
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14),
                      strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
                      strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
                      legend.text = element_text(size = 14),   # legend item labels
                      legend.title = element_text(size = 16)   # legend title
                    )
                }
                
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
                  labs(title = "", x = "Time", y = "Count")+
                  # Increase axis text size
                  theme(
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
                    strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
                    legend.text = element_text(size = 14),   # legend item labels
                    legend.title = element_text(size = 16)   # legend title
                  )
              }
              
              ggsave(
                filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,
                                  "/",prefix,
                                  "param_",param_id,
                                  "_i_opt_",i_opt,
                                  "_overwrite_",overwrite_in,
                                  "_sterile_", sterile_in,
                                  "_tregs_on_",tregs_on_in,
                                  "_tregs_rnd_",trnd_in,
                                  "_",variables[1],"_20.png"),
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
  }
  
  if(skip_id==0){
    # Get unique opt_ind values
    opt_vals = sort(unique(col_avg_keep$opt_ind))
    
    if(skip_opt_0==1){
      # Create color palette:
      n_other = length(opt_vals) 
      viridis_cols = viridis::viridis(n_other)
      
      # Get one extra color and skip the first (lightest) one
      # blues_cols = RColorBrewer::brewer.pal(n_other + 1, "Blues")[-1]
      
      blues_func = colorRampPalette(RColorBrewer::brewer.pal(n_other + 1, "Blues"))
      blues_cols = blues_func(n_other + 1)[-1]
      
      custom_colors = c(setNames(blues_cols, as.character(opt_vals[opt_vals != 0])))
    }else{
      # Create color palette: red for 0, viridis for the rest
      n_other = length(opt_vals) - 1  # number of non-zero values
      viridis_cols = viridis::viridis(n_other)
      
      # Get one extra color and skip the first (lightest) one
      # blues_cols = RColorBrewer::brewer.pal(n_other + 1, "Blues")[-1]
      
      blues_func = colorRampPalette(RColorBrewer::brewer.pal(n_other + 1, "Blues"))
      blues_cols = blues_func(n_other + 1)[-1]
      
      
      custom_colors = c("0" = "red", setNames(blues_cols, as.character(opt_vals[opt_vals != 0])))
    }
    
    # Extract the opt_ind=0, overwrite=0 data
    opt0_ow0 = col_avg_keep %>% filter(opt_ind == 0, overwrite == 0)
    
    # Create a copy with overwrite=1
    opt0_ow1 = opt0_ow0 %>% mutate(overwrite = 1)
    
    # Bind it to your original data
    col_avg_keep = bind_rows(col_avg_keep, opt0_ow1)
    
    p_opt = ggplot(col_avg_keep, aes(x = ros_level, y = avg_pct, color = factor(opt_ind))) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_y_continuous(labels = scales::percent) +
      labs(
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
      facet_wrap(~ overwrite, labeller = label_both, nrow = 1)
    
    ggsave(
      filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,
                        "/A_param_",param_id,
                        "_i_opt_ALL",
                        "_sterile_", sterile_in,
                        "_tregs_on_",tregs_on_in,
                        "_tregs_rnd_",trnd_in,"_20.png"),
      plot = p_opt,
      width = 16,
      height = 8,
      dpi = 300,
      bg='white'
    )
    
    all_nz_inds = unique(col_avg_keep$opt_ind)
    all_nz_inds = all_nz_inds[all_nz_inds>0]
    
    for (i in all_nz_inds){
      
      col_avg_keep_temp = col_avg_keep %>% dplyr::filter(opt_ind %in% c(0,i))
      
      p_opt = ggplot(col_avg_keep_temp, aes(x = ros_level, y = avg_pct, color = factor(opt_ind))) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_y_continuous(labels = scales::percent) +
        labs(
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
        facet_wrap(~ overwrite, labeller = label_both, nrow = 1)
      
      ggsave(
        filename = paste0("/Users/burcutepekule/Desktop/",save_images_path_data,
                          "/A_param_",param_id,
                          "_i_opt_ALL",
                          "_sterile_", sterile_in,
                          "_tregs_on_",tregs_on_in,
                          "_tregs_rnd_",trnd_in,"_20_oi_",i,".png"),
        plot = p_opt,
        width = 16,
        height = 8,
        dpi = 300,
        bg='white'
      )
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
