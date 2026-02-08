
sample_with_efficiency_beta <- function(tau, efficiency, max_concentration = 1000) {
  if (efficiency == 1) return(tau)  # Exact value when perfect efficiency
  if (efficiency == 0) return(runif(1))  # Pure uniform when no efficiency
  
  # Concentration increases exponentially with efficiency
  # This ensures very tight distribution near efficiency = 1
  concentration <- max_concentration * efficiency^3  # Cubic makes it tighter faster
  
  alpha <- 1 + concentration * tau
  beta <- 1 + concentration * (1 - tau)
  
  rbeta(1, shape1 = alpha, shape2 = beta)
}
# Visualize
library(ggplot2)
tau <- 0.7
efficiencies <- c(0, 0.25, 0.5, 0.75, 0.99)

results <- lapply(efficiencies, function(eff) {
  data.frame(
    value = replicate(1000, sample_with_efficiency_beta(tau, eff)),
    efficiency = eff
  )
})

df <- do.call(rbind, results)
df$efficiency <- factor(df$efficiency)

ggplot(df, aes(x = value, fill = efficiency)) +
  geom_density(alpha = 0.7, position = "identity") +
  geom_vline(xintercept = tau, linetype = "dashed", color = "red") +
  facet_wrap(~efficiency, ncol = 1) +
  labs(title = "Distribution as efficiency varies",
       subtitle = paste("tau =", tau))
