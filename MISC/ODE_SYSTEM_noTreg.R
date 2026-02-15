ode_system = function(t, state, params) {
  with(as.list(c(state, params)), {
    
    # Commensals, Lumen (placeholder — fill in later)
    dCL = 0
    
    # Pathogens, Lumen
    dPL = PL*r_L_P*(1-PL/K_L_P) -
      PL*tau_A*(E/E_max)*(1-E/E_max) -
      PL*gamma_B_P*E/E_max
    
    # Pathogens, Lamina Propria
    dPLP = PLP*r_LP_P*(1-PLP/K_LP_P) +
      PL*gamma_B_P*E/E_max -
      PLP*R*r_ros_M - 
      PLP*(1/M_max)*(M_0*alpha_0+M_1*alpha_1+M_2*alpha_2)

    # Commensals, Lamina Propria
    dCLP = CLP*r_LP_C*(1-CLP/K_LP_C) +
      CL*gamma_B_C*E/E_max -
      CLP*R*r_ros_M -
      CLP*(1/M_max)*(M_0*alpha_0+M_1*alpha_1+M_2*alpha_2)
    
    # Epithelial Damage
    dE = PL*tau_P_A*(1-E/E_max) +
      PLP*tau_P_B*(1-E/E_max) + 
      R*r_ros_E*(1-E/E_max) -
      kappa_E*E/E_max
    
    # Danger signals
    dD = beta_D*E/E_max + 
      beta_M*(PLP/K_LP_P + CLP/K_LP_C) +
      beta_A*PLP/K_LP_P -
      tau_D*D
    
    # ROS
    dR = beta_1*M_1/M_max - tau_R*R
    
    # M0 macrophages (placeholder — fill in later)
    dM0 = rec_rate*D*(1-(M_0+M_1+M_2)/M_max) - 
      tau_1_a*M_0/M_max*D + 
      tau_1_d*M_1/M_max*(1/(1+D)) - 
      layoff_rate*(1/(1+D))*M_0/M_max

    # M1 macrophages
    dM1 = tau_1_a*M_0/M_max*D - 
      tau_1_d*M_1/M_max*(1/(1+D))
    
    # M2 macrophages (placeholder — fill in later)
    dM2 = 0
    
    # Order must match state vector: PL, CL, PLP, CLP, E, D, R, M_0, M_1, M_2
    list(c(dPL, dCL, dPLP, dCLP, dE, dD, dR, dM0, dM1, dM2))
  })
}