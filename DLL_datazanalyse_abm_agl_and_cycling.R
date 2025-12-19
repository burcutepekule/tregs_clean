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

acf_peak_val_th   = 0.5
acf_regularity_th = 0.03

score_type     = 'pathogen'
data_suffix    = '_ros_vs_ctrl' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# --- I/O
df_params       = read_csv("./lhs_parameters_della.csv", show_col_types = FALSE)
df_results_keep = readRDS(paste0("./data_cpp_read_abm", data_suffix, ".rds"))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
if(data_suffix == '_ros_vs_ctrl'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else{
  keep_param_id = reps_df %>% dplyr::filter(Freq==100) %>% dplyr::pull(Var1) # 100 = 10 reps per scenario, 5 scenarios x 2 times recording for epithelial and pathogen scores 
}
df_results = df_results_keep %>% dplyr::filter(param_set_id %in% keep_param_id)
print(c(length(unique(df_results_keep$param_set_id)),length(unique(df_results$param_set_id))))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500 # used to be 4500, just for simulation purposes to save time
 
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
num_params = length(unique(df_results$param_set_id))

unique(as.data.frame(table(df_results$param_set_id))[2]) # 40, sanity check!

df_results_keep = df_results
df_results      = df_results %>% dplyr::filter(control==0)
df_results      = df_results %>% 
  dplyr::mutate(
    oscillating = acf_peak_val > acf_peak_val_th & acf_regularity < acf_regularity_th
  )

majority_vote = df_results %>%
  dplyr::group_by(score_type, param_set_id) %>%
  dplyr::summarise(
    n_true = sum(oscillating),
    n_total = n(),
    prop_true = mean(oscillating),
    majority_oscillating = sum(oscillating) > n()/2,
    .groups = "drop"
  )

label_summary = majority_vote %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(
    label_either = any(majority_oscillating), # either epithelial or pathogen
    # label_either = all(majority_oscillating), # both
    oscillating_source = case_when(
      all(majority_oscillating) ~ "both",
      majority_oscillating[score_type == "epithelium"] ~ "epithelium",
      majority_oscillating[score_type == "pathogen"] ~ "pathogen",
      TRUE ~ "none"
    ),
    .groups = "drop"
  )

oscillator_df = label_summary %>%
  dplyr::mutate(
    oscillator_type = case_when(
      label_either ~ "oscillating",
      TRUE ~ "not oscillating"
    ),
    oscillator_type = factor(oscillator_type, levels = c("not oscillating","oscillating"))
  )
oscillator_df = oscillator_df[c('param_set_id','oscillator_type','oscillating_source')]

ss_df      = df_results_keep[c('param_set_id','ss_start','score_type')]
ss_df_mean = ss_df %>%
  dplyr::group_by(param_set_id, score_type) %>%
  dplyr::summarise(
    mean_ss_start = mean(ss_start),
    .groups = "drop"  # removes all grouping
  )

oscillating_merged = merge(ss_df_mean, oscillator_df, by=c('param_set_id'))
parameters_merged  = merge(oscillating_merged, df_params, by=c('param_set_id'))

## =============================================================================

summary_df_param_set_osc = parameters_merged %>% dplyr::filter(oscillator_type=='oscillating')
saveRDS(summary_df_param_set_osc, 'df_agl_cycling_osc.rds')

summary_df_param_set_nosc = parameters_merged %>% dplyr::filter(oscillator_type=='not oscillating')
saveRDS(summary_df_param_set_nosc, 'df_agl_cycling_nosc.rds')

###############################################################################
mult=0
source('~/Dropbox/tregs_clean/MISC/LOAD_PARAM_VECTOR_ROS.R')

df_params_merged  = parameters_merged %>% dplyr::mutate(osc=ifelse(oscillator_type=='oscillating',1,0))

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

library(effsize)
effect_sizes = df_plot_params %>%
  dplyr::group_by(parameter) %>%
  summarise(
    cohens_d = cohen.d(value[osc == "oscillating"], 
                       value[osc == "not oscillating"])$estimate,
    .groups = "drop"
  )

# Create label with effect size category
effect_sizes = effect_sizes %>%
  mutate(
    d_abs = abs(cohens_d),
    effect_label = case_when(
      d_abs < 0.2 ~ paste0("d=", round(cohens_d, 2), " (negligible)"),
      d_abs < 0.5 ~ paste0("d=", round(cohens_d, 2), " (small)"),
      d_abs < 0.8 ~ paste0("d=", round(cohens_d, 2), " (medium)"),
      TRUE ~ paste0("d=", round(cohens_d, 2), " (large)")
    )
  )

# Get y position for each facet (top of each panel)
label_positions = df_plot_params %>%
  group_by(parameter) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  left_join(effect_sizes, by = "parameter")

p_params = ggplot(df_plot_params, aes(x = osc, y = value, fill = osc)) +
  geom_violin(alpha = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~parameter, scales = "free_y", ncol = 4) + # ordered alphabetically (easier to compare among inj types.)
  scale_fill_manual(
    name = "Behaviour",  # or whatever you want
    values = c(
    "not oscillating" = "gray80",
    "oscillating" = "pink"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.10))) +  # Add this line
  geom_text(
    data = label_positions,
    aes(x = 1.5, y = y_pos, label = effect_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  # geom_signif(
  #   comparisons = list(c("not oscillating", "oscillating")),
  #   test = "t.test",
  #   test.args = list(var.equal = FALSE),   # Welch t-test
  #   map_signif_level = TRUE,
  #   textsize = 5,
  #   step_increase = 0.1,
  #   tip_length = 0.02
  # )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  labs(title = paste0("Parameters Distinguishing Oscillatory Behaviour. Number of parameter sets: ",num_params),
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

