library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS('./data_cpp_read_abm.rds')
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
keep_param_id = reps_df %>% dplyr::filter(Freq==200) %>% dplyr::pull(Var1) # 200 = 10 reps per scenario, 10 scenarios x 2 times recording for pathogen scores 
df_results    = df_results_keep %>% filter(param_set_id %in% keep_param_id)
unique(table(df_results$param_set_id)) #66, sanity check
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id)) # >99%!
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold # TRUE, sanity check
unique(table(df_results$param_set_id)) #600, sanity check
length(unique(df_results$param_set_id))

df_comparisons = distinct(df_results %>% dplyr::select(
  param_set_id, injury_type,
  # Mean scores
  mean_ctrl_sterile_e, mean_ctrl_pathogen_e,
  mean_tregs_off_sterile_e, mean_tregs_off_pathogen_e,
  mean_tregs_on_sterile_e, mean_tregs_on_pathogen_e,
  mean_tregs_rnd_sterile_e, mean_tregs_rnd_pathogen_e,
  mean_macspec_sterile_e, mean_macspec_pathogen_e,
  mean_ctrl_sterile_p, mean_ctrl_pathogen_p,
  mean_tregs_off_sterile_p, mean_tregs_off_pathogen_p,
  mean_tregs_on_sterile_p, mean_tregs_on_pathogen_p,
  mean_tregs_rnd_sterile_p, mean_tregs_rnd_pathogen_p,
  mean_macspec_sterile_p, mean_macspec_pathogen_p,
  # SD scores
  sd_ctrl_sterile_e, sd_ctrl_pathogen_e,
  sd_tregs_off_sterile_e, sd_tregs_off_pathogen_e,
  sd_tregs_on_sterile_e, sd_tregs_on_pathogen_e,
  sd_tregs_rnd_sterile_e, sd_tregs_rnd_pathogen_e,
  sd_macspec_sterile_e, sd_macspec_pathogen_e, 
  sd_ctrl_sterile_p, sd_ctrl_pathogen_p,
  sd_tregs_off_sterile_p, sd_tregs_off_pathogen_p,
  sd_tregs_on_sterile_p, sd_tregs_on_pathogen_p,
  sd_tregs_rnd_sterile_p, sd_tregs_rnd_pathogen_p,
  sd_macspec_sterile_p, sd_macspec_pathogen_p
))
df_comparisons_keep = df_comparisons

# ============= FILTER BASED ON CONTROL
df_comparisons_ctrl_test = df_comparisons_keep %>%
  dplyr::mutate(diff_ctrl_test = mean_ctrl_pathogen_p-mean_tregs_off_pathogen_p)

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test[c('param_set_id','diff_ctrl_test')]
df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple, distinct(df_results[c('param_set_id','d_ctrl_vs_test_pathogen_p')]), by='param_set_id')
ids_matter_df = df_comparisons_ctrl_test_simple %>% dplyr::filter(abs(d_ctrl_vs_test_pathogen_p)>=jsd_th & diff_ctrl_test>tol_in)
ids_matter    = unique(ids_matter_df %>% dplyr::pull(param_set_id))
length(ids_matter)
df_comparisons = df_comparisons %>% dplyr::filter(param_set_id %in% ids_matter)

# ============= FILTER BASED ON CONTROL

df_results_pick_p = df_results %>% dplyr::filter(injury_type=='pathogenic' & d_macspec_vs_tregs_on_pathogen_e>0.3)
df_results_pick_s = df_results %>% dplyr::filter(injury_type=='sterile' & d_macspec_vs_tregs_on_sterile_e>0.3)

mean(df_results_pick_p$mean_tregs_on_pathogen_e)
mean(df_results_pick_p$mean_macspec_pathogen_e)

mean(df_results_pick_s$mean_tregs_on_sterile_e)
mean(df_results_pick_s$mean_macspec_sterile_e)


df_comparisons = df_comparisons %>% dplyr::mutate(diff_tregs_on_minus_off = ifelse(injury_type=='pathogenic', 
                                                                                 mean_tregs_on_pathogen_e-mean_macspec_pathogen_e,
                                                                                 mean_tregs_on_sterile_e-mean_macspec_sterile_e))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_tregs_on_minus_off)

df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id',
                                                             'd_macspec_vs_tregs_on_pathogen_e',
                                                             'd_macspec_vs_tregs_on_sterile_e')]), by='param_set_id')

df_comparisons_plot            = df_comparisons

# --- PATHOGENIC ---
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                        ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_macspec_vs_tregs_on_pathogen_e)>jsd_th,
                                                                                                               tregs_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-d_macspec_vs_tregs_on_sterile_e)
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)=='d_macspec_vs_tregs_on_pathogen_e')]='cohens_d'

# --- STERILE ---
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                  ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_macspec_vs_tregs_on_sterile_e)>jsd_th,
                                                                                                         tregs_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-d_macspec_vs_tregs_on_pathogen_e)
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)=='d_macspec_vs_tregs_on_sterile_e')]='cohens_d'

# ============CONFLICT?=========================================================
tregs_better_sterile = df_comparisons_plot_sterile %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_sterile  = df_comparisons_plot_sterile %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)
tregs_better_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::filter(tregs_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_pathogenic  = df_comparisons_plot_pathogenic %>% dplyr::filter(tregs_better_cohens==-1) %>% dplyr::pull(param_set_id)

s_better_p_better = intersect(tregs_better_sterile, tregs_better_pathogenic)
s_better_p_worse  = intersect(tregs_better_sterile, tregs_worse_pathogenic)
s_worse_p_better  = intersect(tregs_worse_sterile, tregs_better_pathogenic)
s_worse_p_worse   = intersect(tregs_worse_sterile, tregs_worse_pathogenic)
# ==============================================================================

df_comparisons_plot = rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)
saveRDS(df_comparisons_plot,'./df_comparisons_plot_macspec.rds')

cohens_th = jsd_th # example threshold for x
e_th      = tol_in # example threshold for y

source('./MISC/REGIONS_macspec.R')

ggsave(
  filename = paste0("./ABM_JSD_MACSPEC.png"),
  plot = p,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)

# print(p)

print(length(unique(dfp$param_set_id)))
dim(df_comparisons_sense)[1] #0! -> ALL MAKES SENSE # IF (sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd), it's insignificant, meaning d_tregs_rnd_vs_nrnd_pathogen<0.5 (or d_tregs_rnd_vs_nrnd_sterile<0.5)!

dfp_p = dfp %>% dplyr::filter(injury_type=='pathogenic') %>% dplyr::select(param_set_id, tregs_better_cohens)
dfp_s = dfp %>% dplyr::filter(injury_type=='sterile') %>% dplyr::select(param_set_id, tregs_better_cohens)
colnames(dfp_p)[2]='tregs_better_cohens_p'
colnames(dfp_s)[2]='tregs_better_cohens_s'
dfp_sp = merge(dfp_p, dfp_s, by='param_set_id')
dfp_sp = dfp_sp %>% dplyr::mutate(treg_better_cohens_both = ifelse(tregs_better_cohens_p+tregs_better_cohens_s==2, 1, 0))

100*length(which(dfp_sp$tregs_better_cohens_p==1))/dim(dfp_sp)[1] # % better for pathogenic
100*length(which(dfp_sp$tregs_better_cohens_p==-1))/dim(dfp_sp)[1] # % worse for pathogenic

100*length(which(dfp_sp$tregs_better_cohens_s==1))/dim(dfp_sp)[1] # % better for sterile
100*length(which(dfp_sp$tregs_better_cohens_s==-1))/dim(dfp_sp)[1] # % worse for sterile

100*sum(dfp_sp$treg_better_cohens_both)/dim(dfp_sp)[1] # % better for both

## check sds?

df_comparisons_plot_nz = df_comparisons_plot %>% dplyr::filter(tregs_better_cohens!=0)

# check conflict?
print(c(s_better_p_worse, s_worse_p_better))
