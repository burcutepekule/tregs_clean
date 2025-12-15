rm(list=ls())
jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = 25*25*0.05
M1_M2_diff     = 1
filter_control = 1
labels_on      = 0
score_type     = 'epithelial' # or 'pathogenic' or 'both'
data_suffix    = '' #
inj_type       = 'sterile'

analysis_pick  = 2 #tregs on vs off
source('./MISC/DECIDE_ANALYSIS.R')
source('./MISC/FILTER_REGIONS.R')
source('./MISC/PLOT_REGIONS.R')
df_comparisons_plot_sc_2 = df_comparisons_plot
round(table(df_comparisons_plot$diff_better_cohens)/dim(df_comparisons_plot)[1],3)
tregs_better_param_ids_2 = unique(df_comparisons_plot_sc_2%>%dplyr::filter(diff_better_cohens==1 & injury_type==inj_type)%>%dplyr::pull(param_set_id))

analysis_pick  = 5 #macspec1 vs tregs on
source('./MISC/DECIDE_ANALYSIS.R')
source('./MISC/FILTER_REGIONS.R')
source('./MISC/PLOT_REGIONS.R')
df_comparisons_plot_sc_5 = df_comparisons_plot
round(table(df_comparisons_plot$diff_better_cohens)/dim(df_comparisons_plot)[1],3)
tregs_better_param_ids_5 = unique(df_comparisons_plot_sc_5%>%dplyr::filter(diff_better_cohens==-1 & injury_type==inj_type)%>%dplyr::pull(param_set_id))

analysis_pick  = 8 #macspec2 vs tregs on
source('./MISC/DECIDE_ANALYSIS.R')
source('./MISC/FILTER_REGIONS.R')
source('./MISC/PLOT_REGIONS.R')
df_comparisons_plot_sc_8 = df_comparisons_plot
round(table(df_comparisons_plot$diff_better_cohens)/dim(df_comparisons_plot)[1],3)
tregs_better_param_ids_8 = unique(df_comparisons_plot_sc_8%>%dplyr::filter(diff_better_cohens==-1 & injury_type==inj_type)%>%dplyr::pull(param_set_id))

intersect(tregs_better_param_ids_2, tregs_better_param_ids_5)
intersect(tregs_better_param_ids_8, tregs_better_param_ids_5)
intersect(tregs_better_param_ids_8, tregs_better_param_ids_2)

setdiff(tregs_better_param_ids_2, tregs_better_param_ids_5)
setdiff(tregs_better_param_ids_5, tregs_better_param_ids_2)

setdiff(tregs_better_param_ids_8, tregs_better_param_ids_5)
setdiff(tregs_better_param_ids_5, tregs_better_param_ids_8)

setdiff(tregs_better_param_ids_8, tregs_better_param_ids_2)
setdiff(tregs_better_param_ids_2, tregs_better_param_ids_8)
