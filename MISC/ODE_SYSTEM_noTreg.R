ode_system = function(t, state, params) {
  with(as.list(c(state, params)), {
    
    # Commensals, Lumen (placeholder — fill in later)
    dCL = 0
    
    # Pathogens, Lumen
    dPL = PL*r_L_P*(1 - PL / K_L_P) -
      PL*f_A(tau_A, E, E_max) -
      PL*f_B(gamma_B_P, E)
    
    # Pathogens, Lamina Propria
    dPLP = PLP*r_LP_P*(1 - PLP / K_LP_P) +
      # PL*f_B(gamma_B_P, E)*(1 - PLP / K_LP_P) -
      PL*f_B(gamma_B_P, E) -
      PLP*f_R(R, theta_M, c_M)*r_ros -
      r_shrink_P*f_shrink(PLP, K_LP_P, c_shrink_P) -
      PLP*(f_E(M_0, alpha_0) + f_E(M_1, alpha_1) + f_E(M_2, alpha_2))
    
    # Commensals, Lamina Propria
    dCLP = CLP*r_LP_C*(1 - CLP / K_LP_C) +
      # CL*f_B(gamma_B_C, E)*(1 - CLP / K_LP_C) -
      CL*f_B(gamma_B_C, E) -
      CLP*f_R(R, theta_M, c_M)*r_ros -
      r_shrink_C*f_shrink(CLP, K_LP_C, c_shrink_C) -
      CLP*(f_E(M_0, alpha_0) + f_E(M_1, alpha_1) + f_E(M_2, alpha_2))
    
    # Epithelial Damage
    dE = PL*tau_P_A*(E_max-E) +
      PLP*tau_P_B*(E_max-E) + f_R(R, theta_E, c_E)*(E_max-E) -
      (E > 0)*kappa_E
    
    # Danger signals
    dD = beta_D*E + beta_M*(PLP + CLP) +
      beta_A*PLP -
      tau_D*D
    
    # ROS
    dR = beta_1*M_1 - tau_R*R
    
    # M0 macrophages (placeholder — fill in later)
    dM0 = rec_rate*D*(M_max-(M_0+M_1)) - 
      tau_1_a*M_0*f_1(D,theta_1,c_1) + 
      tau_1_d*M_1*f_1_rev(D,theta_1,c_1) - 
      layoff_rate*(1/(1+D))*M_0
    
    # M1 macrophages
    dM1 = tau_1_a*M_0*f_1(D,theta_1,c_1) - 
      tau_1_d*M_1*f_1_rev(D,theta_1,c_1)
    
    # M2 macrophages (placeholder — fill in later)
    dM2 = 0
    
    # Order must match state vector: PL, CL, PLP, CLP, E, D, R, M_0, M_1, M_2
    list(c(dPL, dCL, dPLP, dCLP, dE, dD, dR, dM0, dM1, dM2))
  })
}