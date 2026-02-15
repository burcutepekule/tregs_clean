rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(deSolve)
library(stringr)
library(zoo)
library(philentropy)
library(purrr)
library(data.table)
source('./MISC/ODE_SYSTEM_noTreg.R')

# ros_vals   = 1e-2*seq(0,1,0.2)
ros_vals   = c(0, 1, 2, 5, 8, 10, 25, 50, 80)
pat_vals   = c(20, 26, 26.5, 27, 28, 29, 30, 35)
t_max      = 1000
tstep      = 0.1
out_df_all = c()

for (r_ros_in in ros_vals){
  for (r_L_P_in in pat_vals){
    print(c(r_ros_in, r_L_P_in))
    
    params = list(
      
      # Pathogen & Commensal growth (lamina propria)
      # For now, assume no growth is happening in the lamina propria
      r_LP_P    = 0.0,    # P^{LP} logistic growth rate
      r_LP_C    = 0.0,    # C^{LP} logistic growth rate
      K_LP_P    = 1,    # P^{LP} carrying capacity
      K_LP_C    = 1,    # C^{LP} carrying capacity
      
      
      # Pathogen growth (lumen)
      r_L_P     = r_L_P_in,    # P^L logistic growth rate
      K_L_P     = 100,    # P^L carrying capacity
      K_L_C     = 10,    # C^L carrying capacity

      # AMP killing
      tau_A     = 120,    # AMP killing parameter

      # Barrier breach
      gamma_B_P = 0.1,    # barrier breach rate, pathogens
      gamma_B_C = 0.05,   # barrier breach rate, commensals

      # ROS killing of microbes / epithelium
      r_ros_M = 1, # microbes
      r_ros_E = 0.1, # epithelium

      # Engulfment
      alpha_0   = 0,      # engulfment rate by M_0 - assume 0
      alpha_2   = 0,      # engulfment rate by M_2 - assume 0
      alpha_1   = 0.01,    # engulfment rate by M_1

      # Epithelial damage
      tau_P_A   = 0,  # apical damage rate from P^L
      tau_P_B   = 15, # basolateral damage rate from P^{LP}

      # Danger signals
      beta_D    = 1, # DAMP production from epithelial damage
      beta_M    = 0.1, # DAMP production from microbes touching basolateral side
      beta_A    = 0.1, # PAMP production from pathogens in LP
      tau_D     = 0.2, # danger signal decay rate

      # ROS dynamics
      beta_1    = r_ros_in, # ROS production by M_1
      tau_R     = 1,   # ROS decay rate

      # Macrophage dynamics
      theta_1   = 50,    # M1 activation threshold for D
      c_1       = 0.1,  # M1 activation Hill coefficient
      tau_1_a   = 10,   # M1 activation rate
      tau_1_d   = 10,   # M1 deactivation rate

      # Maximum epithelial damage possible
      E_max   = 100, # max epithelial damage
      M_max   = 100, # max number of macrophages
      rec_rate= 10,
      layoff_rate = 1,

      kappa_E     = 1 # recovery rate of epithelium
    )
    
    # --- Initial conditions ---
    state = c(
      PL  = 0.01*params$K_L_P, # start at CC
      CL  = params$K_L_C, 
      PLP = 0,
      CLP = 0,
      E   = 10,
      D   = 0,
      R   = 0,
      M_0 = params$M_max*0.75,
      M_1 = 0,
      M_2 = 0
    )
    
    # --- Solve ---
    times  = seq(0, t_max, by = tstep)
    # out    = as.data.frame(ode(y = state, times = times, func = ode_system, parms = params))
    # Define the extinction threshold
    extinction_threshold = 1e-3  # or whatever makes biological sense
    
    # Root function: triggers when any population crosses threshold
    rootfun = function(t, state, parameters) {
      with(as.list(state), {
        # Return values that become zero at threshold crossings
        c(E - extinction_threshold, 
          PL - extinction_threshold,  
          PLP - extinction_threshold,
          CLP - extinction_threshold,
          M_0 - extinction_threshold,
          M_1 - extinction_threshold)  
      })
    }
    
    # Event function: what to do when root is found
    eventfun = function(t, state, parameters) {
      if (state["E"] < extinction_threshold) state["E"] = 0
      if (state["PL"] < extinction_threshold) state["PL"] = 0
      if (state["PLP"] < extinction_threshold) state["PLP"] = 0
      if (state["CLP"] < extinction_threshold) state["CLP"] = 0
      if (state["M_0"] < extinction_threshold) state["M_0"] = 0
      if (state["M_1"] < extinction_threshold) state["M_1"] = 0
      return(state)
    }
    
    # Then in your ode call:
    out = ode(y = state, times = times, func = ode_system, parms = params,
              rootfun = rootfun, events = list(func = eventfun, root = TRUE))
    
    out_df = as.data.frame(out) 
    out_df$ros_level = r_ros_in
    out_df$pat_level = r_L_P_in
    out_df_all = rbind(out_df_all, out_df)
    
  }
}

variables_2_plot = c('E','PLP','PL','D','R','M_1','M_0')

path_img = './ode_out'
dir.create(path_img)

for (p_ind  in 1:length(variables_2_plot)){
  variables = variables_2_plot[p_ind][[1]]
  variables_name  = paste(variables, collapse = "_")
  
  data_long = out_df_all %>%
    dplyr::select(time, ros_level, pat_level, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p=ggplot(data_long, aes(x = time, y = value))+
    # Lines with agent colors
    geom_line() +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
    theme_minimal() +
    # labs(title = title_opt, x = "Time", y = "Count")
    labs(title = "", x = "Time", y = "Count")+
    # ylim(0, 100) + 
    # Increase axis text size
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text.x = element_text(size = 10),  # top facet labels (ros_level)
      strip.text.y = element_text(size = 10),   # side facet labels (pat_level)
      legend.text = element_text(size = 14),   # legend item labels
      legend.title = element_text(size = 16)   # legend title
    )
  
  ggsave(
    filename = paste0(path_img,"/ode_",variables_name,".png"),
    plot = p,
    width = 18,
    height = 10,
    # width = 12,
    # height = 6,
    dpi = 300,
    bg='white'
  )
}
