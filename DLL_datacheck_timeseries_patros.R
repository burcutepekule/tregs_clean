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

# Get the directory of the current script
get_script_dir = function() {
  # Try rstudioapi first (works in RStudio)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  # Try commandArgs (works with Rscript)
  args = commandArgs(trailingOnly = FALSE)
  file_arg = grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # Try sys.frame (works with source())
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(sys.frame(1)$ofile))
  }
  # Fallback to current working directory
  return(getwd())
}

setwd(get_script_dir())

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Define the ros and pat value ranges
ros_vals = c(0, 0.25, seq(0.5,12,0.5))
pat_vals = c(1:15)
# ros_vals = c(3)
# pat_vals = c(8)
tregs_on_in_vec = c(0)
sterile_in      = 0

# ros_vals = c(0)
# pat_vals = 2
# tregs_on_in_vec = c(0,1)

# ros_vals    = 0
# pat_vals    = 3
# tregs_on_in = 0

# Define both paths
path  = "/Users/burcutepekule/Desktop/sim_abm/"

# ===== reread to include local results 

# param_id_vec = readRDS('evo_selected_triangular_pattern.rds')
param_id_vec = 60800
# TF_matricies = readRDS('control_matrices_all_triangular_patterns.rds')
rep_ind_vec  = 0:4
alpha_plot   = 2/length(rep_ind_vec)

# variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
variables_2_plot = list("epithelial_score","pathogen",c("phagocyte_M1","phagocyte_M2"),c("phagocyte_M1M2","phagocyte_M0"),c("treg_resting", "treg_active"),c("P_M1","P_M2","P_M0"))
# variables_2_plot = list("epithelial_score","pathogen")

background_on    = c(1,1,rep(0,length(variables_2_plot)-2))

# variables_2_plot = list(c("phagocyte_M1","phagocyte_M2","phagocyte_M0"))
# background_on    = 0

for(tregs_on_in in tregs_on_in_vec){
  for(param_id in param_id_vec){
    
    # TF_matrix = TF_matricies[[as.character(param_id)]]
    source('~/Dropbox/tregs_clean/DLL_dataread_reread_patros.R')
    TF_matrix = control_matrix
    # Convert TF_matrix to a data frame for plotting
    TF_df = as.data.frame.table(TF_matrix, stringsAsFactors = FALSE)
    colnames(TF_df) = c("pat_level", "ros_level", "controlled")
    
    # Extract numeric values from row/column names
    TF_df$pat_level = as.numeric(gsub("pat", "", TF_df$pat_level))
    TF_df$ros_level = as.numeric(gsub("ros", "", TF_df$ros_level))
    
    # # Check if output files already exist
    # file1 = paste0("/Users/burcutepekule/Desktop/timeseries_tri/level_", level_in, "_", param_id, "_epithelial_score.png")
    # file2 = paste0("/Users/burcutepekule/Desktop/timeseries_tri/level_", level_in, "_", param_id, "_pathogen.png")
    # 
    # if (file.exists(file1) && file.exists(file2)) {
    #   message("Skipping param_id ", param_id, " - files already exist")
    #   next
    # }
    
    # Dynamically read all RDS files for this parameter set
    # In the loop, check both locations:
    results_list = list()
    for (ros in ros_vals) {
      for (pat in pat_vals) {
        var_name = paste0("results_", ros, "_", pat)
        file_name = paste0('longitudinal_df_param_set_id_', param_id, '_sterile_',sterile_in,'_macspec_0_tregs_',tregs_on_in,'_ros_level_', ros, '_pat_level_', pat, '_trnd_0.rds')
        
        # Try main path first, then local
        file_path = paste0(path, file_name)
        
        if (file.exists(file_path)) {
          results_list[[var_name]] = readRDS(file_path)
        } else {
          warning(paste("File not found in either location:", file_name))
        }
      }
    }
    
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
        geom_hline(yintercept = 150*0.75,
                   linetype = "solid",
                   color = "gray",
                   linewidth = 0.5) +
        # Lines with agent colors
        geom_line(aes(color = variable, group = rep_id), 
                  alpha = alpha_plot, linewidth = 1) +
        scale_color_manual(values = agent_colors, name = "Agent") +
        facet_grid(pat_level ~ ros_level, labeller = label_both) +
        theme_minimal() +
        scale_y_log10(limits = c(NA, 150)) +
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
          labs(title = "", x = "Time", y = "Count")
      }
      
      ggsave(
        filename = paste0("/Users/burcutepekule/Desktop/timeseries_tri/sterile_",sterile_in,"_tregs_on_",tregs_on_in,"_",param_id,"_",variables[1],".png"),
        plot = p,
        width = 24,
        height = 16,
        dpi = 300,
        bg='white'
      )
      
    }
  } 
}

