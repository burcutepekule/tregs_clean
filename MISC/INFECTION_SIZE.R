# Exponential decay from 1 to alpha
inf_load = function(t, alpha_inf, rate_inf) {
  alpha_inf+(1-alpha_inf)*exp(-rate_inf*t)
}

# # Quick check
# f(0)    # Should be 1
# f(Inf)  # Should be 0.1

# Visualize
t = seq(0, 2500, length.out = 1000)

alpha_inf = 0.05
rate_inf  = 0.001

plot(t, inf_load(t, alpha_inf, rate_inf), type = "l", ylim = c(0, 1),
     xlab = "t", ylab = "f(t)",
     main = "Exponential decay from 1 to α")
abline(h = alpha_inf, lty = 2, col = "red")