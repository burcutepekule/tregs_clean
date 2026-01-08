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
path  = "/Users/burcutepekule/Desktop/sim_opt_random/"

# Define the ros and pat value ranges
# ros_vals        = c(3, 5, 10, 20, 50) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
ros_vals        = c(3) # 0 is control - max(ros_level) x max(add_ROS) = 2 x 0.5 = 1 (anyway capped at 1 so makes sense)
pat_vals        = c(10)
iter_vals       = seq(1,100,1)
rep_vals        = seq(1,5,1)
n2_vals         = seq(1,400,1)
tregs_on_in     = 1
sterile_in      = 0
param_id        = 50002

results_list_long  = list()
results_list_theta = list()

path = '/Users/burcutepekule/Desktop/sim_opt_random'

# File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
# Dynamically create file lists for all ros/pat combinations
for (ros in ros_vals) {
  for (pat in pat_vals) {
    for (n2 in n2_vals) {
      for (iter in iter_vals) {
        for (rep in rep_vals) {
          print(paste0("results_", ros, "_", pat, "_", n2, "_", iter, "_", rep))
          file_path_theta = paste0(path,"/theta_param_set_id_",param_id,
                                     "_sterile_0_macspec_0_tregs_1_ros_level_", ros, 
                                     "_pat_level_", pat, 
                                     "_trnd_0_n1_400_n2_",n2,
                                     "_iter_",iter,
                                     "_rep_",rep,".rds")
          
          if(file.exists(file_path_theta)){
            var_name = paste0("results_", ros, "_", pat, "_", n2, "_", iter, "_", rep)
            results_list_theta[[var_name]] = readRDS(file_path_theta)
          }
        }
      }
    }
  }
}

# Combine all results
results_theta = do.call(rbind, results_list_theta)

param_names = c("diffusion_speed_SAMPs",
                 "add_SAMPs",
                 "SAMPs_decay",
                 "treg_discrimination_efficiency",
                 "activation_threshold_SAMPs")

results_theta_df = results_theta %>%
  as.data.frame() %>%
  tibble::rownames_to_column("id") %>%
  separate(id, into = c("prefix", "ros", "pat", "n2", "iter", "rep"), 
           sep = "_", convert = TRUE) %>%
  select(-prefix)

# Rename V columns using the param_names vector
names(results_theta_df)[6:10] = param_names
results_theta_df = results_theta_df %>% dplyr::mutate(unique_id = paste0('ros_',ros,'_pat_',pat,'_n2_',n2,'_iter_',iter))

opt_theta = readRDS(paste0('./df_opt_rnd_75_',param_id,'_use.rds'))

matching_rows = results_theta_df %>%
  semi_join(opt_theta, by = c("diffusion_speed_SAMPs", "add_SAMPs", "SAMPs_decay", 
                              "treg_discrimination_efficiency", "activation_threshold_SAMPs"))

### NOW GO BACK AND READ THESE LONGIDUTINAL FILES PRECISELY
for (ros in unique(matching_rows$ros)) {
  for (pat in unique(matching_rows$pat)) {
    for (n2 in unique(matching_rows$n2)) {
      for (iter in unique(matching_rows$iter)) {
        for (rep in unique(matching_rows$rep)) {
          print(paste0("results_", ros, "_", pat, "_", n2, "_", iter, "_", rep))
          file_path_long = paste0(path,"/longitudinal_df_param_set_id_",param_id,
                                  "_sterile_0_macspec_0_tregs_1_ros_level_", ros,
                                  "_pat_level_", pat,
                                  "_trnd_0_n1_400_n2_",n2,
                                  "_iter_",iter,
                                  "_rep_",rep,".rds")
          
          if(file.exists(file_path_long)){
            var_name = paste0("results_", ros, "_", pat, "_", n2, "_", iter, "_", rep)
            results_list_long[[var_name]] = readRDS(file_path_long)
          }
        }
      }
    }
  }
}
results_long  = do.call(rbind, results_list_long)

results_long_df = results_long %>%
  as.data.frame() %>%
  tibble::rownames_to_column("id") %>%
  separate(id, into = c("prefix", "ros", "pat", "n2", "iter", "rep", "t2"),
           sep = "[_.]", convert = TRUE) %>%
  select(-prefix)
# Check if they match
all(results_long_df$t == results_long_df$t2) #TRUE
# If they do, you can drop t2
results_long_df = results_long_df %>% select(-t2)
results_long_df = results_long_df %>% dplyr::mutate(unique_id = paste0('ros_',ros,'_pat_',pat,'_n2_',n2,'_iter_',iter))

# variables = c("epithelial_score")
variables = c("epithelial_score","phagocyte_M1", "phagocyte_M2","commensal","pathogen")

alpha_plot = 1
max_level_injury = 10*25

### go through matches
id_vec = unique(matching_rows$unique_id)

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
    labs(title = title_opt, x = "Time", y = "Count")
  print(p)
  browser()
}




