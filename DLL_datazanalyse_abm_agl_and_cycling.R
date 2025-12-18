rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)
library(ggsignif)
library(ggpubr)

score_type     = 'pathogen'
data_suffix    = '_ros_vs_ctrl' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
if(data_suffix == '_ros_vs_ctrl'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else{
  keep_param_id = reps_df %>% dplyr::filter(Freq==100) %>% dplyr::pull(Var1) # 100 = 10 reps per scenario, 5 scenarios x 2 times recording for epithelial and pathogen scores 
}
df_results = df_results_keep %>% filter(param_set_id %in% keep_param_id)
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500 # used to be 4500, just for simulation purposes to save time
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id)) # >99%!
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold # TRUE, sanity check

merged_df_keep = merge(df_results, df_params, by='param_set_id')
merged_df_keep = merged_df_keep%>%dplyr::filter(control==0 & score_type=='epithelium')

ss_df          = merged_df_keep[c('param_set_id','ss_start')]
ss_df_mean     = ss_df%>%
  group_by(param_set_id) %>%
  summarise(
    mean_ss_start  = mean(ss_start)
  )
merged_df_keep = distinct(merged_df_keep[c('param_set_id','replicate_id','active_age_limit','acf_peak_val','acf_regularity')])

# # Strong, regular oscillator
# good_oscillator = acf_peak_val > 0.5 & acf_regularity < 0.05

summary_df_param_set = merged_df_keep %>%
  group_by(param_set_id) %>%
  summarise(
    mean_acf_regularity  = mean(acf_regularity),
    mean_acf_peak_val    = mean(acf_peak_val)
  )

summary_df_param_set = summary_df_param_set %>% 
  dplyr::mutate(
    best_oscillator = mean_acf_peak_val > 0.5 & mean_acf_regularity < 0.01
  )

summary_df_param_set = summary_df_param_set %>%
  dplyr::mutate(
    oscillator_type = case_when(
      best_oscillator ~ "oscillating",
      TRUE ~ "not oscillating"
    ),
    oscillator_type = factor(oscillator_type, levels = c("not oscillating","oscillating"))
  )

summary_df_param_set = merge(summary_df_param_set, distinct(merged_df_keep[c('param_set_id','active_age_limit')]))

p_osc = ggplot(summary_df_param_set, aes(x = oscillator_type, y = active_age_limit, fill = oscillator_type)) +
  geom_violin(alpha = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("not oscillating" = "grey70","oscillating" = "tomato"), guide = "none") +
  geom_signif(
    comparisons = list(c( "not oscillating","oscillating")),
    test = "t.test",
    test.args = list(var.equal = FALSE),   # Welch t-test
    map_signif_level = TRUE,
    textsize = 5,
    step_increase = 0.1,
    tip_length = 0.02
  )+
  labs(x = "Oscillator Type", y = "Active Age Limit") +
  theme_minimal()

ggsave(
  filename = paste0("./agl_cycling/agl_cycling_dist.png"),
  plot = p_osc,
  width = 10,
  height = 8,
  dpi = 300,
  bg='white'
)


summary_df_param_set = merge(summary_df_param_set, ss_df_mean, by='param_set_id')

summary_df_param_set_osc = summary_df_param_set %>% dplyr::filter(oscillator_type=='oscillating')
saveRDS(summary_df_param_set_osc, 'df_agl_cycling_osc.rds')

summary_df_param_set_nosc = summary_df_param_set %>% dplyr::filter(oscillator_type=='not oscillating')
saveRDS(summary_df_param_set_nosc, 'df_agl_cycling_nosc.rds')

###############################################################################
mult=0
source('~/Dropbox/tregs_clean/MISC/LOAD_PARAM_VECTOR_ROS.R')

df_params_merged  = merge(df_params, summary_df_param_set[c('param_set_id','oscillator_type')], by='param_set_id')
df_params_merged  = df_params_merged %>% dplyr::filter(param_set_id %in% summary_df_param_set$param_set_id)
df_params_merged  = df_params_merged %>% dplyr::mutate(osc=ifelse(oscillator_type=='oscillating',1,0))

df_plot_params = df_params_merged %>%
  dplyr::select(osc, all_of(param_names)) %>%
  pivot_longer(cols = -osc, names_to = "parameter", values_to = "value") %>%
  mutate(
    param_order = match(parameter, param_names),
    parameter_labeled = paste0(sprintf("%02d", param_order), ". ", parameter),
    parameter_labeled = factor(parameter_labeled, 
                               levels = paste0(sprintf("%02d", 1:length(param_names)), ". ", param_names))
  )


df_plot_params = df_plot_params %>% dplyr::mutate(osc=ifelse(osc==0,'not oscillating','oscillating'))

p_params = ggplot(df_plot_params, aes(x = osc, y = value, fill = osc)) +
  geom_violin(alpha = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~parameter, scales = "free_y", ncol = 4) + # ordered alphabetically (easier to compare among inj types.)
  scale_fill_manual(values = c(
    "not oscillating" = "gray80",
    "oscillating" = "pink"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.10))) +  # Add this line
  geom_signif(
    comparisons = list(c("not oscillating", "oscillating")),
    test = "t.test",
    test.args = list(var.equal = FALSE),   # Welch t-test
    map_signif_level = TRUE,
    textsize = 5,
    step_increase = 0.1,
    tip_length = 0.02
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  ) +
  labs(title = "Top Parameters Distinguishing Effect Regions",
       x = "", y = "Parameter Value")

# print(p_params)

ggsave(
  filename = paste0("./agl_cycling/agl_cycling_dist_all_params.png"),
  plot = p_params,
  width = 18,
  height = 20,
  dpi = 300,
  bg='white'
)


# ggplot(summary_df_param_set, aes(x = active_age_limit, y = mean_acf_peak_val)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   stat_cor(method = "pearson") +  # requires ggpubr
#   theme_minimal()
# 
# ggplot(summary_df_param_set, aes(x = active_age_limit, y = mean_acf_regularity)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, col='red') +
#   stat_cor(method = "pearson") +  # requires ggpubr
#   theme_minimal()
# 
# ggplot(summary_df_param_set, aes(x = good_oscillator, y = active_age_limit, fill = good_oscillator)) +
#   geom_boxplot() +
#   scale_x_discrete(labels = c("FALSE" = "No", "TRUE" = "Yes")) +
#   scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "pink"), guide = "none") +
#   labs(x = "Good Oscillator", y = "Active Age Limit") +
#   theme_minimal()
# 
# length(unique(df_results$param_set_id))
# 
# summary_df_param_set_23206 = summary_df_param_set %>% dplyr::filter(param_set_id==23206)
