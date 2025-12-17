rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS_CPP.R")
source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# eff = 0.073
plot_data_stats_keep = c()
for (eff in seq(0,1,0.1)){
  
  # Create data frame for plotting
  plot_data = data.frame(
    real_ratio = rep(seq(0, 1, 0.1), 1000),
    perceived_ratio = NA
  )
  
  for(i in 1:nrow(plot_data)) {
    rat_com_pat_real = plot_data$real_ratio[i]
    
    # =========================================================================
    # Uniform (Option A)
    alpha_A = 1
    beta_A = 1
    # Centered on truth (Option B)
    concentration = 50  # controls how tight at eff=1
    eps = 0.001
    rat_adj = eps+rat_com_pat_real*(1-2*eps)
    alpha_B = rat_adj*concentration
    beta_B  = (1-rat_adj)*concentration
    # Blend between A and B
    alpha = (1 - eff) * alpha_A + eff * alpha_B
    beta  = (1 - eff) * beta_A  + eff * beta_B
    plot_data$perceived_ratio[i] = sample_rbeta(alpha, beta)
    
    # =========================================================================
    precision_treg = 10*eff  # = 1 at eff=0.1
    alpha = (1-eff)*1 + eff*(rat_com_pat_real*precision_treg)
    beta  = (1-eff)*1 + eff*((1-rat_com_pat_real)*precision_treg)
    plot_data$perceived_ratio[i] = sample_rbeta(alpha, beta)
    
    # =========================================================================
    precision_treg = 10*(exp(5*eff))
    alpha = (1-eff)*1 + eff*(rat_com_pat_real*precision_treg)
    beta  = (1-eff)*1 + eff*((1-rat_com_pat_real)*precision_treg)
    plot_data$perceived_ratio[i] = sample_rbeta(alpha, beta)
    
    
    # =========================================================================
    plot_data$perceived_ratio[i] = eff*rat_com_pat_real+(1-eff)*runif(1)
    
  }
  
  if(eff<0.5 & eff>0){
    color_eff = 'pink'
  }else if(eff==0){
    color_eff = 'gray40'
  }else if(eff>=0.5 & eff<1){
    color_eff = 'lightblue'
  }else{
    color_eff = 'green'
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
  
  plot_data_stats = plot_data %>% 
    dplyr::group_by(real_ratio) %>% 
    dplyr::summarise(group_min = min(perceived_ratio),
                     group_mean = mean(perceived_ratio),
                     group_max = max(perceived_ratio)) 
  
  plot_data_stats$eff=eff
  
  plot_data_stats_keep=rbind(plot_data_stats_keep, plot_data_stats)
  
  ggsave(
    filename = paste0("./precision_",100*eff,".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300,
    bg='white'
  )
}


plot_data_stats_keep = plot_data_stats_keep %>% dplyr::mutate(diff_ratio=group_mean-real_ratio)
plot_data_stats_keep = plot_data_stats_keep %>% dplyr::mutate(width_ratio=group_max-group_min)


plot(plot_data_stats_keep$eff, plot_data_stats_keep$width_ratio)

plot(plot_data_stats_keep$eff, plot_data_stats_keep$diff_ratio)

# print(p)
# 
# plot_data_pick = plot_data %>% dplyr::filter(real_ratio==0.5)
# hist(plot_data_pick$perceived_ratio, 50)
