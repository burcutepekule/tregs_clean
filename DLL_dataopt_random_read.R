library(dplyr)
library(tidyr)
library(zoo)
param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

loop_over_all   = 92000 # c(30000, 81250, 45500, 92000, 62500)
df_opt_rnd      = read.table(paste0("/Users/burcutepekule/Dropbox/tregs_clean/merged_",loop_over_all,".txt"), 
                             header = FALSE, 
                             col.names = c("param_set_id" ,"rep", "pat_level", "ros_level", "overwrite_in",
                                           param_names, "min_e", "min_p","mean_e","mean_p",
                                           "pct_above_threshold_e","pct_above_threshold_p"))

summary_df = df_opt_rnd %>%
  group_by(param_set_id, pat_level, ros_level, diffusion_speed_SAMPs, overwrite_in,
           add_SAMPs, SAMPs_decay, treg_discrimination_efficiency, activation_threshold_SAMPs
           ) %>%
  summarise(
    n_reps = n(),
    mean_pct_above_threshold_e = mean(pct_above_threshold_e),
    mean_pct_above_threshold_p = mean(pct_above_threshold_p),
    .groups = "drop"
  )

summary_df_keep = summary_df

summary_df_095 = summary_df %>% dplyr::filter(mean_pct_above_threshold_e>0.95 & mean_pct_above_threshold_p > 0.95 & n_reps>=3)
df_opt_rnd_95  = summary_df_095[param_names]
saveRDS(summary_df_095,paste0('./summary_df_095_',loop_over_all,'_use.rds'))
saveRDS(df_opt_rnd_95,paste0('./df_opt_rnd_95_',loop_over_all,'_use.rds'))

summary_df_080 = summary_df %>% dplyr::filter(mean_pct_above_threshold_e>0.80 & mean_pct_above_threshold_p > 0.80 & n_reps>=3)
df_opt_rnd_80  = summary_df_080[param_names]
saveRDS(df_opt_rnd_80,paste0('./df_opt_rnd_80_',loop_over_all,'_use.rds'))

summary_df_075 = summary_df %>% dplyr::filter(mean_pct_above_threshold_e>0.75 & mean_pct_above_threshold_p > 0.75 & n_reps>=3)
df_opt_rnd_75  = summary_df_075[param_names]
saveRDS(df_opt_rnd_75,paste0('./df_opt_rnd_75_',loop_over_all,'_use.rds'))

# df_opt_rnd_95  = df_opt_rnd %>% dplyr::filter(pct_above_threshold_e>0.95 & pct_above_threshold_p>0.95 & min_e==150 & min_p==0) 
# df_opt_rnd_95  = df_opt_rnd_95[param_names]
# print(df_opt_rnd_95)

# saveRDS(df_opt_rnd_95,paste0('./df_opt_rnd_95_',loop_over_all,'_2.rds'))

# pick_idxs     = c(15, 26, 31, 32, 34, 41, 24, 38, 47, 22, 43, 28, 37, 49)
# df_opt_rnd_95 = df_opt_rnd_95[pick_idxs,]
# rownames(df_opt_rnd_95) = 1:dim(df_opt_rnd_95)[1]
# 
# saveRDS(df_opt_rnd_95,paste0('./df_opt_rnd_95_',loop_over_all,'_use.rds'))
print(summary_df_095)

