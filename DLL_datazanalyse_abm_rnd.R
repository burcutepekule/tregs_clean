rm(list=ls())
jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = 5*tol_in_e
M1_M2_diff     = 1
filter_control = 0
labels_on      = 0
score_type     = 'epithelial' # or 'pathogenic' or 'both'
# score_type     = 'pathogen' # or 'pathogen' or 'both'
# score_type     = 'both'
# data_suffix    = '_10' # empty for 100 reps, _10 for 10 reps
data_suffix    = '_100' # empty for 100 reps, _10 for 10 reps

inj_type= 'sterile'
# inj_type= 'pathogenic'
# inj_type= 'pooled'

source('./MISC/FILTER_REGIONS.R') #df_comparisons
df_comparisons_keep = df_comparisons

condition_subt_from = 'tregs_on'
condition_subt      = 'tregs_off'
jensen_distance     = 'tregs_on_vs_off'
df_comparisons_in   = df_comparisons_keep
source('./MISC/FILTER_FOR_SUBTRACT.R')
df_comparisons_plot_1 = df_comparisons_plot
tregs_better_when_on  = df_comparisons_plot_1 %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)

condition_subt_from = 'tregs_on'
condition_subt      = 'tregs_rnd'
jensen_distance     = 'tregs_on_vs_rnd'
df_comparisons_in   = df_comparisons_keep
source('./MISC/FILTER_FOR_SUBTRACT.R')
df_comparisons_plot_2 = df_comparisons_plot
tregs_better_when_not_random  = df_comparisons_plot_2 %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_better_when_random      = df_comparisons_plot_2 %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)


tregs_better_on_but_get_worse_when_randomized = intersect(tregs_better_when_on, tregs_better_when_not_random)

setdiff(tregs_better_on_but_get_worse_when_randomized, tregs_better_when_not_random)

setdiff(tregs_better_when_on, tregs_better_when_not_random)
setdiff(tregs_better_when_not_random, tregs_better_when_on)

### NOT ALL PARAM SETS LEAD TO WORSE OUTCOMES WHEN TREGS ARE BENEFICIAL BUT RANDOMIZED
### BUT ALL PARAM SETS THAT ARE BETTER WHEN NOT RANDOMIZED ARE ALSO BETTER WHEN ON!


### Any param sets where tregs are better when random and also beneficial when on?
intersect(tregs_better_when_on, tregs_better_when_random) ##ZERO! YES!

### Do some tregs get







