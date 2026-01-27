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
source("./DLL_dataopt_random_read.R")

path = "/Users/burcutepekule/Desktop/sim_opt_random/"

param_id  = 92000
# param_id   = 30000
params_opt = readRDS(paste0('./summary_df_035_',param_id,'_use.rds'))
# ros_vals  = unique(params_opt$ros_level)
# pat_vals  = unique(params_opt$pat_level)
ros_vals  = unique(params_opt$ros_level)
pat_vals  = unique(params_opt$pat_level)

iter_vals = seq(1, 100, 1)
rep_vals  = seq(1, 5, 1)
n1        = 400
n2_vals   = seq(1, n1, 1)

# 1. Load the optimal theta FIRST
# opt_theta = readRDS(paste0('./df_opt_rnd_75_', param_id, '_use.rds'))
opt_theta = readRDS(paste0('./df_opt_rnd_35_', param_id, '_use.rds'))
print(paste0("opt_theta read"))

# 2. Build file index and read ONLY theta files (smaller), then filter
# Read theta files in chunks or one by one, checking for matches immediately

for(overwrite_in in c(0,1)){
  results_list_long = list()
  results_list_theta = list()
  
  print(paste0("loop starting"))
  for (ros in ros_vals) {
    for (pat in pat_vals) {
      for (n2 in n2_vals) {
        for (iter in iter_vals) {
          for (rep in rep_vals) {
            
            file_path_theta = paste0(path, "/theta_param_set_id_", param_id,
                                     "_overwrite_",overwrite_in,
                                     "_sterile_0_macspec_0_tregs_1_ros_level_", ros,
                                     "_pat_level_", pat,
                                     "_trnd_0_n1_", n1,
                                     "_n2_", n2,
                                     "_iter_", iter,
                                     "_rep_", rep, ".rds")
            
            if (file.exists(file_path_theta)) {
              theta_data = readRDS(file_path_theta)
              theta_df = as.data.frame(t(theta_data))
              names(theta_df) = c("diffusion_speed_SAMPs", "add_SAMPs", "SAMPs_decay",
                                  "treg_discrimination_efficiency", "activation_threshold_SAMPs")
              print(theta_df)
              
              # Check if this theta matches opt_theta IMMEDIATELY
              is_match = nrow(semi_join(theta_df, opt_theta, 
                                        by = c("diffusion_speed_SAMPs", "add_SAMPs", "SAMPs_decay",
                                               "treg_discrimination_efficiency", "activation_threshold_SAMPs"))) > 0
              
              if (is_match) {
                print(paste0("MATCH FOUND: ", ros, "_", pat, "_", n2, "_", iter, "_", rep))
                
                var_name = paste0("results_", ros, "_", pat, "_", n2, "_", iter, "_", rep)
                
                # Store theta
                theta_df$ros = ros
                theta_df$pat = pat
                theta_df$n2 = n2
                theta_df$iter = iter
                theta_df$rep = rep
                theta_df$unique_id = paste0('ros_', ros, '_pat_', pat, '_n2_', n2, '_iter_', iter)
                results_list_theta[[var_name]] = theta_df
                
                # Read longitudinal file ONLY for matches
                file_path_long = paste0(path, "/longitudinal_df_param_set_id_", param_id,
                                        "_overwrite_",overwrite_in,
                                        "_sterile_0_macspec_0_tregs_1_ros_level_", ros,
                                        "_pat_level_", pat,
                                        "_trnd_0_n1_", n1,
                                        "_n2_", n2,
                                        "_iter_", iter,
                                        "_rep_", rep, ".rds")
                
                if (file.exists(file_path_long)) {
                  results_list_long[[var_name]] = readRDS(file_path_long)
                }
              }else{
                print(paste0("MATCH NOT FOUND: ", ros, "_", pat, "_", n2, "_", iter, "_", rep))
              }
            }
          }
        }
      }
    }
  }
  
  # Combine only the matching results
  matching_rows = do.call(rbind, results_list_theta)
  results_long = do.call(rbind, results_list_long)
  
  # Continue with processing
  # results_long_df = results_long %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column("id") %>%
  #   separate(id, into = c("prefix", "ros", "pat", "n2", "iter", "rep", "t2"),
  #            sep = "[_.]", convert = TRUE) %>%
  #   select(-prefix)
  
  results_long_df = results_long %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    separate(id, into = c("prefix", "ros", "pat", "n2", "iter", "rep_t2"),
             sep = "_", convert = FALSE) %>%  # Don't auto-convert yet
    separate(rep_t2, into = c("rep", "t2"), 
             sep = "\\.", convert = TRUE) %>%  # Split rep.t2 on the dot
    mutate(
      ros = as.numeric(ros),
      pat = as.numeric(pat),
      n2 = as.numeric(n2),
      iter = as.numeric(iter)
    ) %>%
    select(-prefix)
  
  # Check if they match
  all(results_long_df$t == results_long_df$t2) #TRUE
  # If they do, you can drop t2
  results_long_df = results_long_df %>% select(-t2)
  results_long_df = results_long_df %>% dplyr::mutate(unique_id = paste0('ros_',ros,'_pat_',pat,'_n2_',n2,'_iter_',iter))
  
  # variables = c("epithelial_score")
  variables = c("epithelial_score","phagocyte_M1", "phagocyte_M2","commensal","pathogen")
  
  alpha_plot = 0.5
  max_level_injury = 10*25
  
  ### go through matches
  id_vec = unique(matching_rows$unique_id)
  dir.create(paste0("/Users/burcutepekule/Desktop/opt_",param_id))
  
  for(id in id_vec){
    matching_rows_temp = matching_rows %>% dplyr::filter(unique_id==id)
    title_opt = paste(distinct(matching_rows_temp[param_names]), collapse = "_")
    print(id)
    print(distinct(matching_rows_temp[param_names]))
    results_long_df_use = results_long_df %>% dplyr::filter(unique_id==id)
    
    data_long = results_long_df_use %>%
      dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    p = ggplot(data_long, aes(x = t, y = value))+
      # Lines with agent colors
      geom_line(aes(color = variable, group = interaction(rep_id, variable)), 
                alpha = alpha_plot, linewidth = 1) +
      scale_color_manual(values = agent_colors, name = "Agent") +
      facet_wrap(~ variable, ncol = 2, labeller = label_both, scales = "free_y") +
      theme_minimal() +
      scale_y_log10(limits = c(NA, max_level_injury)) +
      scale_y_log10() +
      labs(title = paste0("ros_",unique(matching_rows_temp$ros),"_pat_",unique(matching_rows_temp$pat)),
           subtitle = title_opt, 
           x = "Time", y = "Count")
    ggsave(
      filename = paste0("/Users/burcutepekule/Desktop/opt_",param_id,"/i_opt_",id,"_ow_",overwrite_in,".png"),
      plot = p,
      width = 24,
      height = 16,
      dpi = 300,
      bg='white'
    )
    
  }
  
  matching_rows_distinct = distinct(matching_rows %>% dplyr::select(-rep))
  matching_rows_distinct[1,]
}

