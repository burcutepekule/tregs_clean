total_injury= sum(epithelium$level_injury)
amp_killing = amp_kill_rate*(total_injury/(grid_size*max_level_injury))*(1-total_injury/(grid_size*max_level_injury)) # both scaled by the level of injury and leftover healthy cells 

pat_lumen = pat_lumen+pat_lumen*r_pat*(1-pat_lumen/pat_lumen_max)
pat_lumen = pat_lumen-pat_lumen*amp_killing

# Pathogens leak through injured epithelium
pathogen_source       = pat_lumen*epithelium$level_injury*p_leak_constant # this is because I don't wanna change the rate_leak_pathogen_injury in ASSIGN_PARAMETERS.R
density_pathogen[1, ] = density_pathogen[1, ]+pathogen_source
pat_lumen             = pat_lumen - sum(pathogen_source)

# Commensals: baseline leak+injury-enhanced leak (commensals are infinite)
commensal_source_baseline = rep(rate_leak_commensal_baseline*0.01)
commensal_source_injury   = epithelium$level_injury*rate_leak_commensal_injury*0.01
density_commensal[1, ]    = density_commensal[1, ]+commensal_source_baseline+commensal_source_injury