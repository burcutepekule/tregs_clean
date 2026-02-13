library(dplyr)
library(tidyr)
library(zoo)
param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")
loop_over_all = c(172)
m2on_in = 0

for (loop_over in loop_over_all){
  # loop_over = 90 # c(92000, 30000, 81250, 45500, 62500)
  min_reps      = 10 # max is 10
  # df_opt_rnd    = read.table(paste0("/Users/burcutepekule/Dropbox/tregs_clean/merged_",loop_over,"_m2on_",m2on_in,".txt"),
  df_opt_rnd    = read.table(paste0("/Users/burcutepekule/Dropbox/tregs_clean/merged_",loop_over,".txt"),
                             header = FALSE,
                             col.names = c("param_set_id" ,"rep", "pat_level", "ros_level", 
                                           "overwrite_in",param_names, "min_e", "min_p",
                                           "mean_e","mean_p","pct_above_threshold_e","pct_above_threshold_p",
                                           "sum_e","sum_p")) ### NOW OPTIMIZE ACCORDING TO THIS!
  
  summary_df = df_opt_rnd %>%
    group_by(param_set_id, pat_level, ros_level, diffusion_speed_SAMPs, overwrite_in,
             add_SAMPs, SAMPs_decay, treg_discrimination_efficiency, activation_threshold_SAMPs
    ) %>%
    summarise(
      n_reps = n(),
      mean_pct_above_threshold_e = mean(pct_above_threshold_e),
      mean_pct_above_threshold_p = mean(pct_above_threshold_p),
      mean_sum_e = mean(sum_e),
      mean_sum_p = mean(sum_p),
      .groups = "drop"
    )
  
  summary_df_keep  = summary_df
  summary_df_10    = summary_df %>% dplyr::filter(n_reps>=min_reps)
  
  summary_df_10 = summary_df_10 %>% dplyr::rowwise() %>% dplyr::mutate(mean_pct_above_threshold_min = min(mean_pct_above_threshold_e,mean_pct_above_threshold_p))

  saveRDS(summary_df_10,paste0('./summary_df_10rep_',loop_over,'_use_m2on_',m2on_in,'.rds'))
  summary_df_10 = summary_df_10[order(summary_df_10$mean_sum_p, decreasing = FALSE),]
  
  summary_df_10 = summary_df_10[order(summary_df_10$pat_level,summary_df_10$mean_pct_above_threshold_min, decreasing = TRUE),]
  
  summary_df_10_only_1 = summary_df_10 %>% dplyr::filter(mean_pct_above_threshold_e>0.8 & mean_pct_above_threshold_p>0.8)
  print(head(summary_df_10[, c('param_set_id','pat_level','ros_level','mean_pct_above_threshold_e', 'mean_pct_above_threshold_p','overwrite_in','n_reps')], 5))
  
  summary_df_10 = summary_df_10[order(summary_df_10$mean_sum_p, decreasing = FALSE),]
}
