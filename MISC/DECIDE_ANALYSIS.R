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
