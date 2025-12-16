rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# eff = 0.073
eff = .5

precision_treg = 10*(exp(5*eff)) #ORIGINAL
precision_treg = 10*eff  # = 1 at eff=0.1
# precision_treg = exp(5*eff)  # ~1.65 at eff=0.1

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

if(eff<0.5 & eff>0){
  color_eff = 'pink'
}else if(eff==0){
  color_eff = 'gray40'
}else{
  color_eff = 'blue'
}

p=ggplot(plot_data, aes(x = factor(real_ratio), y = perceived_ratio)) +
  geom_violin(fill = color_eff, alpha = 0.2, trim = TRUE) +
  geom_boxplot(fill = color_eff, width = 0.2, alpha = 0.8, outlier.shape = NA) +
  labs(
    x = "Real ratio",
    y = "Distribution of ratio perceived by Treg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

plot_data_means = plot_data %>% 
  dplyr::group_by(real_ratio) %>% 
  dplyr::summarise(group_mean = mean(perceived_ratio))

ggsave(
  filename = paste0("./precision_",eff,".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  bg='white'
)

print(p)