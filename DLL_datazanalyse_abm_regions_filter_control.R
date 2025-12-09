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
keep_param_id = reps_df %>% dplyr::filter(Freq==100) %>% dplyr::pull(Var1) # 100 = 10 reps per scenario, 10 scenarios
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
  param_set_id, injury_type, control, macspec_on,
  # Mean scores
  mean_ctrl_sterile, mean_ctrl_pathogen,
  mean_tregs_off_sterile, mean_tregs_off_pathogen,
  mean_tregs_on_sterile, mean_tregs_on_pathogen,
  mean_tregs_rnd_sterile, mean_tregs_rnd_pathogen,
  mean_macspec_sterile, mean_macspec_pathogen,
  # SD scores
  sd_ctrl_sterile, sd_ctrl_pathogen,
  sd_tregs_off_sterile, sd_tregs_off_pathogen,
  sd_tregs_on_sterile, sd_tregs_on_pathogen,
  sd_tregs_rnd_sterile, sd_tregs_rnd_pathogen,
  sd_macspec_sterile, sd_macspec_pathogen
))
df_comparisons_keep = df_comparisons

### ======= HERE - SHOULD BE PATHOGEN LOAD!==========================================
df_comparisons_ctrl_test = df_comparisons_keep %>% dplyr::filter(macspec_on==0) %>% dplyr::mutate(diff_ctrl_test = ifelse(injury_type=='pathogenic', 
                                                                                                                          mean_ctrl_pathogen-mean_tregs_off_pathogen,
                                                                                                                          mean_ctrl_sterile-mean_tregs_off_sterile))
df_comparisons_ctrl_test_simple = distinct(df_comparisons_ctrl_test[c('param_set_id','injury_type','diff_ctrl_test')])

df_comparisons_ctrl_test_simple = merge(df_comparisons_ctrl_test_simple, distinct(df_results[c('param_set_id','d_ctrl_vs_test_pathogen','d_ctrl_vs_test_sterile')]), by='param_set_id')

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% dplyr::mutate(d_ctrl_test = ifelse(injury_type=='pathogenic', 
                                                                                                         d_ctrl_vs_test_pathogen,
                                                                                                         d_ctrl_vs_test_sterile))

df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple[c('param_set_id','injury_type','diff_ctrl_test','d_ctrl_test')]
df_comparisons_ctrl_test_simple = df_comparisons_ctrl_test_simple %>% dplyr::filter(abs(d_ctrl_test)>jsd_th & abs(diff_ctrl_test)>tol_in)

param_ids_twice = df_comparisons_ctrl_test_simple %>%
  count(param_set_id) %>%
  filter(n == 2) %>%
  pull(param_set_id)
### ======= HERE - SHOULD BE PATHOGEN LOAD!==========================================


d_tregs_on_vs_off_pathogen

df_comparisons = df_comparisons %>% dplyr::filter(control==0 & macspec_on==0)
df_comparisons = df_comparisons %>% dplyr::mutate(diff_tregs_on_minus_off = ifelse(injury_type=='pathogenic', 
                                                                                   mean_tregs_on_pathogen-mean_tregs_off_pathogen,
                                                                                   mean_tregs_on_sterile-mean_tregs_off_sterile))

df_comparisons = df_comparisons %>% dplyr::mutate(diff_notrnd_minus_rnd = ifelse(injury_type=='pathogenic', 
                                                                                 mean_tregs_on_pathogen-mean_tregs_rnd_pathogen,
                                                                                 mean_tregs_on_sterile-mean_tregs_rnd_sterile))

df_comparisons = df_comparisons %>% dplyr::mutate(tregs_off_sd = ifelse(injury_type=='pathogenic', sd_tregs_off_pathogen, sd_tregs_off_sterile))
df_comparisons = df_comparisons %>% dplyr::mutate(tregs_on_no_rnd_sd = ifelse(injury_type=='pathogenic', sd_tregs_on_pathogen, sd_tregs_on_sterile))
df_comparisons = df_comparisons %>% dplyr::mutate(tregs_on_rnd_sd = ifelse(injury_type=='pathogenic', sd_tregs_rnd_pathogen, sd_tregs_rnd_sterile))

df_comparisons = df_comparisons %>% dplyr::select(param_set_id, injury_type, diff_tregs_on_minus_off, diff_notrnd_minus_rnd, 
                                                  tregs_off_sd, tregs_on_no_rnd_sd, tregs_on_rnd_sd)

df_comparisons = merge(df_comparisons, distinct(df_results[c('param_set_id',
                                                             'd_tregs_on_vs_off_pathogen',
                                                             'd_tregs_rnd_vs_nrnd_pathogen',
                                                             'd_tregs_on_vs_off_sterile',
                                                             'd_tregs_rnd_vs_nrnd_sterile')]), by='param_set_id')

#----- filter based on cohen's d since you are looking for differences!
df_comparisons_below= df_comparisons %>% dplyr::filter((injury_type=='pathogenic' & abs(d_tregs_on_vs_off_pathogen)<jsd_th)|(injury_type=='sterile' & abs(d_tregs_on_vs_off_sterile)<jsd_th))
# tol_in              = 125*0.15 # this is NOT over epithelial cells (so don't use 25*stg), it's over epithelial score! max diff is 6*25-1*25=125!
# or get tol_in based on the range that is not significant!
if(is.na(tol_in)){
  tol_in = max(max(df_comparisons_below$diff_tregs_on_minus_off), abs(min(df_comparisons_below$diff_tregs_on_minus_off)))
}

# # === CUTOFF FOR VARIANCE? STABLE PARAMETER SETS?
# log10_variance = log10(c(df_comparisons$tregs_off_sd, df_comparisons$tregs_on_no_rnd_sd, df_comparisons$tregs_on_rnd_sd))
# d = density(log10_variance)
# y = d$y
# x = d$x
# # Local maxima (peaks)
# peaks_idx = which(diff(sign(diff(y))) == -2) + 1
# peaks_x   = x[peaks_idx]
# peaks_y   = y[peaks_idx]
# # Local minima (valleys)
# valleys_idx = which(diff(sign(diff(y))) == 2) + 1
# valleys_x   = x[valleys_idx]
# valleys_y   = y[valleys_idx]
# 
# plot(x, y, type="l")
# points(peaks_x, peaks_y, col="red", pch=19, cex=1.5)
# points(valleys_x, valleys_y, col="blue", pch=19, cex=1.5)
# 
# sd_tol_in = 10^valleys_x[2] # take the second, less conservative?

# Your variable (epithelial score) spans 100 units: between 25 to 125 
# For a distribution on that scale:
# SD ≈ 5–10 → very low variation
# SD ≈ 10–20 → moderate variation
# SD ≈ 20–30+ → high variation

#----- cases where randomizing tregs have the opposite effect to turning them on, and significant?

df_comparisons_keep  = df_comparisons
df_comparisons_sense = df_comparisons %>% dplyr::mutate(makes_sense=ifelse(sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd) 
                                                                           & abs(diff_notrnd_minus_rnd)>tol_in
                                                                           & (abs(d_tregs_rnd_vs_nrnd_pathogen)>=jsd_th | abs(d_tregs_rnd_vs_nrnd_sterile)>=jsd_th) ,0,1)) #either for sterile or pathogenic
df_comparisons_sense = df_comparisons_sense %>% dplyr::filter(makes_sense==0) 
dim(df_comparisons_sense)[1] #0! -> ALL MAKES SENSE # IF (sign(diff_tregs_on_minus_off)!=sign(diff_notrnd_minus_rnd), it's insignificant, meaning d_tregs_rnd_vs_nrnd_pathogen<0.5 (or d_tregs_rnd_vs_nrnd_sterile<0.5)!

#----- Find significant cases

df_comparisons_plot            = df_comparisons_keep
# --- PATHOGENIC ---
df_comparisons_plot_pathogenic = df_comparisons_plot %>% dplyr::filter(injury_type=='pathogenic')
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                        ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_tregs_on_vs_off_pathogen)>jsd_th,
                                                                                                               # & tregs_on_no_rnd_sd<sd_tol_in
                                                                                                               # & tregs_off_sd<sd_tol_in,
                                                                                                               tregs_better, 0))

df_comparisons_plot_pathogenic = df_comparisons_plot_pathogenic %>% dplyr::select(-d_tregs_on_vs_off_sterile)
colnames(df_comparisons_plot_pathogenic)[which(colnames(df_comparisons_plot_pathogenic)=='d_tregs_on_vs_off_pathogen')]='cohens_d'

# --- STERILE ---
df_comparisons_plot_sterile = df_comparisons_plot %>% dplyr::filter(injury_type=='sterile')
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better = ifelse(diff_tregs_on_minus_off > tol_in, 1,
                                                                                                  ifelse(diff_tregs_on_minus_off < -1*tol_in,-1,0)))
df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::mutate(tregs_better_cohens = ifelse(abs(d_tregs_on_vs_off_sterile)>jsd_th,
                                                                                                         # & tregs_on_no_rnd_sd<sd_tol_in
                                                                                                         # & tregs_off_sd<sd_tol_in,
                                                                                                         tregs_better, 0))

df_comparisons_plot_sterile = df_comparisons_plot_sterile %>% dplyr::select(-d_tregs_on_vs_off_pathogen)
colnames(df_comparisons_plot_sterile)[which(colnames(df_comparisons_plot_sterile)=='d_tregs_on_vs_off_sterile')]='cohens_d'

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
saveRDS(df_comparisons_plot,'./df_comparisons_plot.rds')

cohens_th = jsd_th # example threshold for x
e_th      = tol_in # example threshold for y

source('./MISC/REGIONS.R')

ggsave(
  filename = paste0("./ABM_JSD_filtered.png"),
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
