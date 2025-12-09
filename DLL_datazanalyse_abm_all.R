rm(list=ls())
sd_tol_in   = Inf 
jsd_th      = 0.3
tol_in      = 125*0.2
M1_M2_diff  = 1
filter_control = 0
source('./DLL_datazanalyse_abm_regions_tregs_on_vs_tregs_off.R')

df_params           = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep     = readRDS('./data_cpp_read_abm.rds')
df_plot             = readRDS('./df_comparisons_plot.rds')
df_plot             = merge(df_plot, df_params, by='param_set_id')
if(M1_M2_diff==1){
  df_plot = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
}
source('./MISC/LOAD_PARAM_VECTOR.R') #M1_M2_diff adjusts params as well

inj_type= 'sterile'
inj_type= 'pathogenic'
# inj_type= 'pooled'

if(inj_type!='pooled'){
  df_plot = df_plot %>% dplyr::filter(injury_type==inj_type)
}
df_lda  = df_plot %>% dplyr::select(all_of(param_names), tregs_better_cohens)
classes = unique(df_lda$tregs_better_cohens)
 
# Choose confidence levels for visualization
level_plus1  = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
level_minus1 = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
violin_on    = 1
if(identical(classes, c(0,1))){
  source('./MISC/PLS_DA_2classes.R')
}else{
  source('./MISC/PLS_DA_3classes.R')
}
