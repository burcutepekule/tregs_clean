rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

rate_leak_pathogen_injury = .1

r_pat         = rate_leak_pathogen_injury
pat_lumen     = 10
pat_lumen_max = pat_lumen
inj_level     = c(rep(0,5),rep(1,15),rep(0,5))

p_leak_constant=(1/(250))
pathogen_source_all = c()
pat_lumen_all = c()

for (t in 1:100){
  pat_lumen           = pat_lumen*(1+r_pat*(1-pat_lumen/pat_lumen_max))
  pathogen_source     = pat_lumen*inj_level*p_leak_constant # this is because I don't wanna change the rate_leak_pathogen_injury in ASSIGN_PARAMETERS.R
  pat_lumen           = pmax(0,pat_lumen - sum(pathogen_source))
  pathogen_source_all = c(pathogen_source_all, sum(pathogen_source))
  pat_lumen_all       = c(pat_lumen_all, pat_lumen)
}

plot(pat_lumen_all)
plot(pathogen_source_all)
