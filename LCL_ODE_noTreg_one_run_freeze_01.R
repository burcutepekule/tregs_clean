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
source('./MISC/ODE_FUNCTIONS.R')
source('./MISC/ODE_SYSTEM_noTreg.R')

t_max   = 5000

params = list(
  # Pathogen growth (lumen)
  r_L_P     = 0.2,    # P^L logistic growth rate
  K_L_P     = 100,    # P^L carrying capacity
  K_L_C     = 10,    # C^L carrying capacity
  
  # Pathogen & Commensal growth (lamina propria)
  # For now, assume no growth is happening in the lamina propria
  r_LP_P    = 0.0,    # P^{LP} logistic growth rate
  r_LP_C    = 0.0,    # C^{LP} logistic growth rate
  K_LP_P    = 1,    # P^{LP} carrying capacity
  K_LP_C    = 1,    # C^{LP} carrying capacity
  c_shrink_P = 0.1, 
  c_shrink_C = 0.1,
  r_shrink_P = 1, #on/off
  r_shrink_C = 1,
  
  # AMP killing
  tau_A     = 1.08,    # AMP killing parameter
  
  # Barrier breach
  gamma_B_P = 0.0001,    # barrier breach rate, pathogens
  gamma_B_C = 0.00005,   # barrier breach rate, commensals
  
  # ROS killing of microbes
  theta_M   = 250, # ROS threshold (microbes)
  c_M       = 0.05, # ROS killing Hill coefficient (microbes)
  r_ros     = 0.01, # killing rate
  
  # Engulfment
  alpha_0   = 0,      # engulfment rate by M_0 - assume 0 
  alpha_2   = 0,      # engulfment rate by M_2 - assume 0 
  alpha_1   = 0.1,    # engulfment rate by M_1
  
  # Epithelial damage
  tau_P_A   = 0.0001,  # apical damage rate from P^L
  tau_P_B   = 0.001, # basolateral damage rate from P^{LP}
  theta_E   = 100,   # ROS damage threshold
  c_E       = 0.1,    # ROS damage Hill coefficient (epithelium)
  kappa_E   = 0.05,    # epithelial recovery rate
  
  # Danger signals
  beta_D    = 0.01, # DAMP production from epithelial damage
  beta_M    = 0.001, # DAMP production from microbes touching basolateral side
  beta_A    = 0.001, # PAMP production from pathogens in LP
  tau_D     = 0.002, # danger signal decay rate
  
  # ROS dynamics
  beta_1    = 0.0,    # ROS production by M_1
  tau_R     = 0.01,   # ROS decay rate
  
  # Macrophage dynamics
  theta_1   = 50,    # M1 activation threshold for D
  c_1       = 0.1,  # M1 activation Hill coefficient
  tau_1_a   = 0.001,   # M1 activation rate
  tau_1_d   = 0.001,   # M1 deactivation rate
  
  # Maximum epithelial damage possible
  E_max   = 100, # max epithelial damage 
  M_max   = 100, # max number of macrophages
  rec_rate= 0.001,
  layoff_rate = 0.001
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
times  = seq(0, t_max, by = 0.1)
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

out_df = out_df[c('time','PL','PLP','CLP','E','M_0','M_1','D','R')]
col_order = setdiff(names(out_df), "time")
out_df %>%
  pivot_longer(-time, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = col_order)) %>%
  ggplot(aes(x = time, y = value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal()