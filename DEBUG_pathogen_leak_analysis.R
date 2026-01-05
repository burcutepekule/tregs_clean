# ============================================================================
# DIAGNOSTIC SCRIPT: Why does rate_leak_pathogen_injury have variable effects?
# ============================================================================
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the simulation results
results_dir = '/scratch/gpfs/CMETCALF/sim_abm'
param_set_id = 60800

# Get all RDS files for this parameter set
files = list.files(results_dir,
                   pattern = paste0("longitudinal_df_param_set_id_", param_set_id, ".*\\.rds$"),
                   full.names = TRUE)

cat("Found", length(files), "simulation files for parameter set", param_set_id, "\n")

# Read all files and combine
all_data = lapply(files, function(f) {
  df = readRDS(f)

  # Extract scenario parameters from filename
  fname = basename(f)
  parts = strsplit(fname, "_")[[1]]

  df$sterile = as.numeric(gsub("sterile", "", parts[grep("sterile", parts)]))
  df$macspec = as.numeric(gsub("macspec", "", parts[grep("macspec", parts)]))
  df$tregs = as.numeric(gsub("tregs", "", parts[grep("tregs", parts)]))
  df$ros_level = as.numeric(gsub("level", "", parts[grep("ros.*level", parts)]))
  df$pat_level = as.numeric(gsub("level", "", parts[grep("pat.*level", parts)]))

  return(df)
}) %>% bind_rows()

# Calculate summary statistics
summary_stats = all_data %>%
  group_by(ros_level, pat_level, rep_id) %>%
  summarize(
    # Epithelial outcomes (last 500 timesteps = steady state)
    mean_epithelial_score_ss = mean(epithelial_score[t >= 1500], na.rm=TRUE),
    sd_epithelial_score_ss = sd(epithelial_score[t >= 1500], na.rm=TRUE),

    # Pathogen outcomes
    mean_pathogen_ss = mean(pathogen[t >= 1500], na.rm=TRUE),
    max_pathogen = max(pathogen, na.rm=TRUE),

    # Immune response metrics
    mean_M1_ss = mean(phagocyte_M1[t >= 1500], na.rm=TRUE),
    mean_ROS_ss = mean(P_ROS[t >= 1500], na.rm=TRUE),

    .groups = 'drop'
  )

# Average across replicates
summary_by_scenario = summary_stats %>%
  group_by(ros_level, pat_level) %>%
  summarize(
    epithelial_score = mean(mean_epithelial_score_ss, na.rm=TRUE),
    epithelial_score_se = sd(mean_epithelial_score_ss, na.rm=TRUE) / sqrt(n()),
    pathogen = mean(mean_pathogen_ss, na.rm=TRUE),
    M1 = mean(mean_M1_ss, na.rm=TRUE),
    ROS = mean(mean_ROS_ss, na.rm=TRUE),
    .groups = 'drop'
  )

cat("\n=== Analysis Results ===\n\n")

# For each ROS level, analyze the effect of pathogen leak rate
for (ros in sort(unique(summary_by_scenario$ros_level))) {
  cat("ROS level:", ros, "\n")

  subset_data = summary_by_scenario %>% filter(ros_level == ros) %>% arrange(pat_level)

  if (nrow(subset_data) > 1) {
    # Calculate correlation between pat_level and epithelial_score
    cor_val = cor(subset_data$pat_level, subset_data$epithelial_score, use="complete.obs")

    # Calculate slope of epithelial_score vs pat_level
    model = lm(epithelial_score ~ pat_level, data=subset_data)
    slope = coef(model)[2]

    cat("  Correlation (pat_level vs epithelial_score):", round(cor_val, 3), "\n")
    cat("  Slope (epithelial_score per unit pat_level):", round(slope, 3), "\n")

    # Show range of epithelial scores
    cat("  Epithelial score range:",
        round(min(subset_data$epithelial_score, na.rm=TRUE), 1), "to",
        round(max(subset_data$epithelial_score, na.rm=TRUE), 1), "\n")

    # Show immune response strength
    cat("  Mean M1 macrophages:", round(mean(subset_data$M1, na.rm=TRUE), 1), "\n")
    cat("  Mean ROS:", round(mean(subset_data$ROS, na.rm=TRUE), 3), "\n\n")
  }
}

# Create visualization
cat("\n=== Creating visualization ===\n")

# Plot 1: Epithelial score vs pathogen leak rate, faceted by ROS level
p1 = ggplot(summary_by_scenario, aes(x=pat_level, y=epithelial_score, color=factor(ros_level))) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  facet_wrap(~ros_level, ncol=3, labeller = label_both) +
  labs(
    title = "Effect of Pathogen Leak Rate on Epithelial Score",
    subtitle = paste("Parameter set", param_set_id),
    x = "Pathogen leak rate (pat_level)",
    y = "Epithelial score (higher = healthier)",
    color = "ROS level"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("DEBUG_epithelial_score_vs_pat_level.png", p1, width=14, height=10, dpi=300)
cat("Saved: DEBUG_epithelial_score_vs_pat_level.png\n")

# Plot 2: Heatmap showing epithelial score
p2 = ggplot(summary_by_scenario, aes(x=pat_level, y=ros_level, fill=epithelial_score)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "yellow", high = "green",
    midpoint = 90, # halfway between min (25) and max (150)
    name = "Epithelial\nScore"
  ) +
  labs(
    title = "Epithelial Score Heatmap",
    subtitle = "ROS level vs Pathogen leak rate",
    x = "Pathogen leak rate (pat_level)",
    y = "ROS level"
  ) +
  theme_bw()

ggsave("DEBUG_epithelial_heatmap.png", p2, width=10, height=8, dpi=300)
cat("Saved: DEBUG_epithelial_heatmap.png\n")

# Plot 3: Scatter showing the interaction
p3 = ggplot(summary_by_scenario, aes(x=pathogen, y=epithelial_score, color=ros_level, size=M1)) +
  geom_point(alpha=0.6) +
  scale_color_viridis_c() +
  labs(
    title = "Pathogen Load vs Epithelial Score",
    subtitle = "Point size = M1 macrophages",
    x = "Mean pathogen count (steady state)",
    y = "Epithelial score",
    color = "ROS level"
  ) +
  theme_bw()

ggsave("DEBUG_pathogen_vs_epithelial.png", p3, width=10, height=8, dpi=300)
cat("Saved: DEBUG_pathogen_vs_epithelial.png\n")

cat("\n=== Analysis complete! ===\n")
cat("\nKey findings:\n")
cat("1. Check if high ROS levels compensate for high pathogen leak rates\n")
cat("2. Look for parameter combinations where immune response saturates\n")
cat("3. Identify threshold effects where system switches from clearance to collapse\n")
