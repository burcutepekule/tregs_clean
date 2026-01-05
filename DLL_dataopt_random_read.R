rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
param_names = c("diffusion_speed_SAMPs",
                "add_SAMPs",
                "SAMPs_decay",
                "treg_discrimination_efficiency",
                "activation_threshold_SAMPs")

loop_over_all   = 5504
df_opt_rnd      = read.table(paste0("./merged_",loop_over_all,".txt"), 
                             header = FALSE, 
                             col.names = c("param_set_id" ,"rep", "pat_level", "ros_level",
                                           param_names, "min_e", "min_p","mean_e","mean_p",
                                           "pct_above_threshold_e","pct_above_threshold_p"))

df_opt_rnd_95  = df_opt_rnd %>% dplyr::filter(pct_above_threshold_e>0.95 & pct_above_threshold_p>0.95 & min_e==150 & min_p==0) 
df_opt_rnd_95  = df_opt_rnd_95[param_names]
print(df_opt_rnd_95)

# saveRDS(df_opt_rnd_95,paste0('./df_opt_rnd_95_',loop_over_all,'_2.rds'))

# pick_idxs     = c(15, 26, 31, 32, 34, 41, 24, 38, 47, 22, 43, 28, 37, 49)
# df_opt_rnd_95 = df_opt_rnd_95[pick_idxs,]
# rownames(df_opt_rnd_95) = 1:dim(df_opt_rnd_95)[1]
# 
# saveRDS(df_opt_rnd_95,paste0('./df_opt_rnd_95_',loop_over_all,'_use.rds'))
