# flog_th = function(x, th, c){
#   out = ifelse(x > th,
#                1/(1+exp(-c*(x -th))), # might need to play around with 0.1
#                0)
#   return(out)
# }
flog_th = function(x, th, c) {
  s = 1 / (1 + exp(-c * (x - th)))
  ifelse(x > th, 2 * (s - 0.5), 0)  # maps [0.5, 1] → [0, 1]
}

flog_th_rev = function(x, th, c) {
  s = 1-1 / (1 + exp(-c * (x - th)))
  ifelse(x < th, 2 * (s - 0.5), 0)  # maps [0.5, 1] → [0, 1]
}

f_A = function(tau_A, E, E_max) {
  out = tau_A*(E/E_max)*(1-E/E_max)
  return(out)
}
f_B = function(gamma_B, E) {
  out = E*gamma_B
  return(out)
}
f_R = function(R, theta_M, c_M) {
  return(flog_th(R, theta_M, c_M))
}
f_E = function(M_i, alpha_i) {
  return(M_i*alpha_i)
}
f_1 = function(D, theta_1, c_1) {
  return(flog_th(D, theta_1, c_1))
}
f_1_rev = function(D, theta_1, c_1) {
  return(flog_th_rev(D, theta_1, c_1))
}
f_shrink = function(PLP, K_LP_P, c){
  ifelse(PLP>K_LP_P, c*(PLP-K_LP_P), 0)
}