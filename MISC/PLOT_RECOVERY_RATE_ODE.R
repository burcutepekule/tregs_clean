k = 100
epith_recovery_chance = 0.05
E_max = 1
E = seq(0,1,0.001)
recovery_rate = epith_recovery_chance*(1-(exp(k*(E-E_max))))

# recovery_rate = epith_recovery_chance / (1e-6 + exp(k*(E-E_max)))
plot(E, recovery_rate)
