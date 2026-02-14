# --- Parameters ---
params = list(
  # Pathogen growth (lumen)
  # r_L_P     = 0.02,    # P^L logistic growth rate
  r_L_P     = r_L_P_in,
  K_L_P     = 10,    # P^L carrying capacity
  K_L_C     = 1,    # C^L carrying capacity
  
  # Pathogen & Commensal growth (lamina propria)
  # For now, assume no growth is happening in the lamina propria
  r_LP_P    = 0.0,    # P^{LP} logistic growth rate
  r_LP_C    = 0.0,    # C^{LP} logistic growth rate
  K_LP_P    = 10,    # P^{LP} carrying capacity
  K_LP_C    = Inf,    # C^{LP} carrying capacity
  
  # AMP killing
  tau_A     = 0.0,    # AMP killing parameter
  
  # Barrier breach
  gamma_B_P = 0.001,    # barrier breach rate, pathogens
  gamma_B_C = 0.0005,   # barrier breach rate, commensals
  
  # ROS killing of microbes
  theta_M   = 250, # ROS threshold (microbes)
  c_M       = 0.05, # ROS killing Hill coefficient (microbes)
  r_ros     = 0.01, # killing rate

  # Engulfment
  alpha_0   = 0,      # engulfment rate by M_0 - assume 0 
  alpha_2   = 0.01,      # engulfment rate by M_2 - assume 0 
  alpha_1   = 0,    # engulfment rate by M_1
  
  # Epithelial damage
  tau_P_A   = 0.01,  # apical damage rate from P^L
  tau_P_B   = 0.01, # basolateral damage rate from P^{LP}
  theta_E   = 100,   # ROS damage threshold
  c_E       = 0.1,    # ROS damage Hill coefficient (epithelium)
  kappa_E   = 0.01,    # epithelial recovery rate
  
  # Danger signals
  beta_D    = 0.01, # DAMP production from epithelial damage
  beta_M    = 0.001, # DAMP production from microbes touching basolateral side
  beta_A    = 0.001, # PAMP production from pathogens in LP
  tau_D     = 0.001, # danger signal decay rate
  
  # ROS dynamics
  beta_1    = r_ros_in,    # ROS production by M_1
  tau_R     = 0.01,   # ROS decay rate
  
  # Macrophage dynamics
  theta_1   = 0.1,    # M1 activation threshold
  c_1       = 0.1,  # M1 activation Hill coefficient
  tau_1_a   = 0.1,   # M1 activation rate
  tau_1_d   = 0.05,   # M1 deactivation rate
  
  # Maximum epithelial damage possible
  E_max   = 100, # max epithelial damage 
  M_max   = 100, # max number of macrophages
  rec_rate= 0.01,
  layoff_rate = 0
)

# --- Initial conditions ---
state = c(
  PL  = params$K_L_P, # start at CC
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