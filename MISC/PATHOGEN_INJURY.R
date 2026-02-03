rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

pat_vec = seq(0,1000,0.1)


k_in  = 0.01
x0_in = 5
f     = x0_in*(1-exp(-k_in*pat_vec))

plot(pat_vec, (f))
lines(pat_vec, ceiling(f), col='purple')
ceiling(f)

# Logistic function parameters (for epithelial injury calculation)
k_in  = 0.45
x0_in = 10
c_in  = 10
lines(logistic_scaled_0_to_5_quantized(pat_vec), col='blue')
logistic_scaled_0_to_5_quantized(0)
logistic_scaled_0_to_5_quantized(1)

# Logistic function parameters (for epithelial injury calculation)
k_in_log  = 0.044
x0_in_log = 50
c_in_log  = 5
plot(logistic_scaled_0_to_5_quantized(pat_vec,k_in_log, x0_in_log, c_in_log), col='red')
lines(logistic_scaled_0_to_5_quantized(pat_vec,k_in_log, x0_in_log, c_in_log), col='red')


