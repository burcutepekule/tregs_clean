rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

eff = .1
precision_treg = 10*(exp(5*eff*eff))

# Create data frame for plotting
plot_data = data.frame(
  real_ratio = rep(seq(0, 1, 0.1), 1000),
  perceived_ratio = NA
)

for(i in 1:nrow(plot_data)) {
  rat_com_pat_real = plot_data$real_ratio[i]
  alpha = (1-eff)*1 + eff*(rat_com_pat_real*precision_treg)
  beta  = (1-eff)*1 + eff*((1-rat_com_pat_real)*precision_treg)
  plot_data$perceived_ratio[i] = sample_rbeta(alpha, beta)
}

ggplot(plot_data, aes(x = factor(real_ratio), y = perceived_ratio)) +
  geom_boxplot(fill = "pink", alpha = 0.7) +
  labs(
    x = "Real ratio",
    y = "Distribution of ratio perceived by Treg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
