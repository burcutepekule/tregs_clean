rm(list=ls())
jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = 25*25*0.5
M1_M2_diff     = 1
filter_control = 0
labels_on      = 0
score_type     = 'epithelial' # or 'pathogenic' or 'both'
# score_type     = 'pathogen' # or 'pathogen' or 'both'
# score_type     = 'both'
# data_suffix    = '_10' # empty for 100 reps, _10 for 10 reps
# data_suffix    = '_100' # empty for 100 reps, _10 for 10 reps
data_suffix    = '' # empty for 100 reps, _10 for 10 reps

inj_type= 'sterile'
inj_type= 'pathogenic'
# inj_type= 'pooled'

source('./MISC/FILTER_REGIONS.R') #df_comparisons
df_comparisons_keep = df_comparisons %>% dplyr::filter(injury_type==inj_type)

condition_subt_from = 'tregs_on'
condition_subt      = 'tregs_off'
jensen_distance     = 'tregs_on_vs_off'
df_comparisons_in   = df_comparisons_keep
source('./MISC/FILTER_FOR_SUBTRACT.R')
df_comparisons_plot_1 = df_comparisons_plot
tregs_better_when_on  = df_comparisons_plot_1 %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_worse_when_on   = df_comparisons_plot_1 %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)

condition_subt_from = 'tregs_on'
condition_subt      = 'tregs_rnd'
jensen_distance     = 'tregs_on_vs_rnd'
df_comparisons_in   = df_comparisons_keep
source('./MISC/FILTER_FOR_SUBTRACT.R')
df_comparisons_plot_2 = df_comparisons_plot
tregs_better_when_not_random  = df_comparisons_plot_2 %>% dplyr::filter(diff_better_cohens==1) %>% dplyr::pull(param_set_id)
tregs_better_when_random      = df_comparisons_plot_2 %>% dplyr::filter(diff_better_cohens==-1) %>% dplyr::pull(param_set_id)
tregs_dm_when_random          = df_comparisons_plot_2 %>% dplyr::filter(diff_better_cohens==0) %>% dplyr::pull(param_set_id)

tregs_worse_when_not_random  = tregs_better_when_random
tregs_worse_when_random      = tregs_better_when_not_random

tregs_better_on_but_get_worse_when_randomized = intersect(tregs_better_when_on, tregs_worse_when_random)

p1=length(tregs_better_on_but_get_worse_when_randomized)/length(tregs_better_when_on)

### 82% of the cases where Tregs have a benefit lose that benefit when randomized!
## What about the rest?

tregs_better_on_and_perform_similar_when_randomized = intersect(tregs_better_when_on, tregs_dm_when_random)
p2=length(tregs_better_on_and_perform_similar_when_randomized)/length(tregs_better_when_on)
## 17% does not get affected by randomization

p1+p2 #[1] 0.9949495

## what is left? benefit INCREASING when randomized 
setdiff(tregs_better_when_on, c(tregs_better_on_but_get_worse_when_randomized,tregs_better_on_and_perform_similar_when_randomized))

## this is tregs_better_on_but_get_better_when_randomized -> # there is only 1!
## can also be calculated tregs_better_on_but_get_better_when_randomized = intersect(tregs_better_when_on, tregs_better_when_random)

########
# tregs_worse_on_but_get_better_when_randomized = intersect(tregs_worse_when_on, tregs_better_when_random)
# # actually tregs_worse_on_but_get_better_when_randomized is not contradictory.
# # what would be contradictory is tregs are better compared to off, and better comprared to non-random. 







