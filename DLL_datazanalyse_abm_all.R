rm(list=ls())
jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = 25*25*0.05
M1_M2_diff     = 1
filter_control = 1
labels_on      = 0
score_type     = 'epithelial' # or 'pathogenic' or 'both'
# score_type     = 'pathogen' # or 'pathogenic' or 'both'

data_suffix    = '' #
# data_suffix    = '_pre_pamps' #

### THIS IS FOR PARAMETERS
inj_type= 'sterile' 
inj_type= 'pathogenic'
inj_type= 'pooled'

analysis_pick  = 2 #2, 5, 8?

if(analysis_pick==1){
  # 1 =========================================
  condition_subt_from = 'ctrl'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'ctrl_vs_tregs_off'
}else if(analysis_pick==2){
  # 2 =========================================
  condition_subt_from = 'tregs_on'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'tregs_on_vs_off'
}else if(analysis_pick==3){
  # 3 =========================================
  condition_subt_from = 'tregs_on'
  condition_subt      = 'tregs_rnd'
  jensen_distance     = 'tregs_on_vs_rnd'
}else if(analysis_pick==4){
  # 4 =========================================
  condition_subt_from = 'macspec1'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'macspec1_vs_tregs_off'
}else if(analysis_pick==5){
  # 5 =========================================
  condition_subt_from = 'macspec1'
  condition_subt      = 'tregs_on'
  jensen_distance     = 'macspec1_vs_tregs_on'
}else if(analysis_pick==6){
  # 6 =========================================
  condition_subt_from = 'macspec1'
  condition_subt      = 'tregs_rnd'
  jensen_distance     = 'macspec1_vs_tregs_rnd'
}else if(analysis_pick==7){
  # 7 =========================================
  condition_subt_from = 'macspec2'
  condition_subt      = 'tregs_off'
  jensen_distance     = 'macspec2_vs_tregs_off'
}else if(analysis_pick==8){
  # 8 =========================================
  condition_subt_from = 'macspec2'
  condition_subt      = 'tregs_on'
  jensen_distance     = 'macspec2_vs_tregs_on'
}else if(analysis_pick==9){
  # 9 =========================================
  condition_subt_from = 'macspec2'
  condition_subt      = 'tregs_rnd'
  jensen_distance     = 'macspec2_vs_tregs_rnd'
}

source('./MISC/FILTER_REGIONS.R')
source('./MISC/PLOT_REGIONS.R') # pathogen scores reversed
# source('./MISC/PLOT_REGIONS_SIMPLE.R')  # pathogen scores NOT reversed, interpret like epithelium


df_params           = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
colnames(df_params)[which(colnames(df_params)=='th_ROS_epith_recover')] = 'th_ROS_epith_injury'
df_results_keep     = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
df_plot             = df_comparisons_plot
df_plot             = merge(df_plot, df_params, by='param_set_id')
if(M1_M2_diff==1){
  df_plot = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
}
source('./MISC/LOAD_PARAM_VECTOR.R') #M1_M2_diff adjusts params as well

if(inj_type!='pooled'){
  df_plot = df_plot %>% dplyr::filter(injury_type==inj_type)
}
df_lda  = df_plot %>% dplyr::select(all_of(param_names), diff_better_cohens)
classes = unique(df_lda$diff_better_cohens)

# Choose confidence levels for visualization
level_plus1  = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
level_minus1 = 0.75 # pick from c(0.50, 0.75, 0.90, 0.95, 0.99)
violin_on    = 1
if(sum(classes)==1){ 
  source('./MISC/PLS_DA_2classes_0p1.R')
}else if(sum(classes)==-1){
  source('./MISC/PLS_DA_2classes_0n1.R')
}else{
  source('./MISC/PLS_DA_3classes.R')
}
