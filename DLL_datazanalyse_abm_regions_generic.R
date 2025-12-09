# =============================================================================
# Generic ABM Region Analysis Script
# =============================================================================
# This script performs pairwise comparisons between different treatment conditions
# and generates region plots with Cohen's d analysis.
#
# Usage:
# Set the configuration parameters below, then source this script.
# =============================================================================

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Define comparison groups
# Available options: "ctrl", "tregs_off", "tregs_on", "tregs_rnd", "macspec"
group1 <- "macspec"      # First group for comparison
group2 <- "tregs_off"    # Second group for comparison

# Difference calculation: group1 - group2 or group2 - group1
# Set to "group1_minus_group2" or "group2_minus_group1"
difference_direction <- "group1_minus_group2"

# Control filtering
filter_control <- 0  # Set to 1 to filter based on control, 0 to skip

# Thresholds (if not defined in environment, set defaults)
if (!exists("jsd_th")) jsd_th <- 0.3    # Cohen's d threshold
if (!exists("tol_in")) tol_in <- 18.75  # Tolerance for difference

# Output file prefix (will generate: <prefix>.rds and ABM_JSD_<prefix>.png)
output_prefix <- paste0(group1, "_vs_", group2)

# Additional analyses flags
check_randomization_logic <- FALSE  # Set TRUE for tregs_on vs tregs_off analysis
show_intermediate_picks   <- FALSE  # Set TRUE to show intermediate filtering results

# =============================================================================
# LOAD LIBRARIES
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# =============================================================================
# DATA LOADING AND INITIAL FILTERING
# =============================================================================

df_params       <- read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep <- readRDS('./data_cpp_read_abm.rds')
length(unique(df_results_keep$param_set_id))

# Filter for complete number of reps
reps_df       <- as.data.frame(table(df_results_keep$param_set_id))
keep_param_id <- reps_df %>%
  dplyr::filter(Freq == 200) %>%
  dplyr::pull(Var1)
df_results <- df_results_keep %>%
  filter(param_set_id %in% keep_param_id)

cat("After filtering for complete reps:\n")
cat("  Unique param sets:", length(unique(df_results$param_set_id)), "\n")

# Filter based on ss_start threshold
ss_start_threshold <- 4500
param_id_all_below <- df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)

cat("After filtering for ss_start < ", ss_start_threshold, ":\n", sep = "")
cat("  Proportion kept:",
    round(length(param_id_all_below) / length(unique(df_results$param_set_id)) * 100, 1),
    "%\n")

df_results <- df_results %>%
  dplyr::filter(param_set_id %in% param_id_all_below)

cat("  Unique param sets:", length(unique(df_results$param_set_id)), "\n")

# =============================================================================
# PREPARE COMPARISON DATA
# =============================================================================

df_comparisons <- distinct(df_results %>% dplyr::select(
  param_set_id, injury_type,
  # Mean scores
  mean_ctrl_sterile_e, mean_ctrl_pathogen_e,
  mean_tregs_off_sterile_e, mean_tregs_off_pathogen_e,
  mean_tregs_on_sterile_e, mean_tregs_on_pathogen_e,
  mean_tregs_rnd_sterile_e, mean_tregs_rnd_pathogen_e,
  mean_macspec_sterile_e, mean_macspec_pathogen_e,
  mean_ctrl_sterile_p, mean_ctrl_pathogen_p,
  mean_tregs_off_sterile_p, mean_tregs_off_pathogen_p,
  mean_tregs_on_sterile_p, mean_tregs_on_pathogen_p,
  mean_tregs_rnd_sterile_p, mean_tregs_rnd_pathogen_p,
  mean_macspec_sterile_p, mean_macspec_pathogen_p
))

df_comparisons_keep <- df_comparisons

# =============================================================================
# FILTER BASED ON CONTROL (OPTIONAL)
# =============================================================================

if (filter_control == 1) {
  cat("\nApplying control filter...\n")

  df_comparisons_ctrl_test <- df_comparisons_keep %>%
    dplyr::mutate(diff_ctrl_test = mean_ctrl_pathogen_p - mean_tregs_off_pathogen_p)

  df_comparisons_ctrl_test_simple <- df_comparisons_ctrl_test[c('param_set_id', 'diff_ctrl_test')]
  df_comparisons_ctrl_test_simple <- merge(
    df_comparisons_ctrl_test_simple,
    distinct(df_results[c('param_set_id', 'd_ctrl_vs_test_pathogen_p')]),
    by = 'param_set_id'
  )

  ids_matter_df <- df_comparisons_ctrl_test_simple %>%
    dplyr::filter(abs(d_ctrl_vs_test_pathogen_p) >= jsd_th & diff_ctrl_test > tol_in)
  ids_matter <- unique(ids_matter_df %>% dplyr::pull(param_set_id))

  cat("  Param sets passing control filter:", length(ids_matter), "\n")

  df_comparisons <- df_comparisons %>%
    dplyr::filter(param_set_id %in% ids_matter)
}

# =============================================================================
# INTERMEDIATE PICKS (OPTIONAL)
# =============================================================================

if (show_intermediate_picks) {
  cat("\nIntermediate picks analysis:\n")

  # Construct Cohen's d column names
  cohens_col_p <- paste0("d_", group1, "_vs_", group2, "_pathogen_e")
  cohens_col_s <- paste0("d_", group1, "_vs_", group2, "_sterile_e")

  # Check if these columns exist, if not try reverse
  if (!cohens_col_p %in% colnames(df_results)) {
    cohens_col_p <- paste0("d_", group2, "_vs_", group1, "_pathogen_e")
    cohens_col_s <- paste0("d_", group2, "_vs_", group1, "_sterile_e")
  }

  if (cohens_col_p %in% colnames(df_results) && cohens_col_s %in% colnames(df_results)) {
    df_results_pick_p <- df_results %>%
      dplyr::filter(injury_type == 'pathogenic' & abs(.data[[cohens_col_p]]) > 0.3)
    df_results_pick_s <- df_results %>%
      dplyr::filter(injury_type == 'sterile' & abs(.data[[cohens_col_s]]) > 0.3)

    cat("  Pathogenic cases with |Cohen's d| > 0.3:", nrow(df_results_pick_p), "\n")
    cat("  Sterile cases with |Cohen's d| > 0.3:", nrow(df_results_pick_s), "\n")
  }
}

# =============================================================================
# COMPUTE DIFFERENCES BETWEEN GROUPS
# =============================================================================

cat("\nComputing differences between", group1, "and", group2, "...\n")

# Construct variable names dynamically
mean_g1_pathogen_e <- paste0("mean_", group1, "_pathogen_e")
mean_g1_sterile_e  <- paste0("mean_", group1, "_sterile_e")
mean_g2_pathogen_e <- paste0("mean_", group2, "_pathogen_e")
mean_g2_sterile_e  <- paste0("mean_", group2, "_sterile_e")

# Calculate difference based on direction
if (difference_direction == "group1_minus_group2") {
  df_comparisons <- df_comparisons %>%
    dplyr::mutate(
      diff_tregs_on_minus_off = ifelse(
        injury_type == 'pathogenic',
        .data[[mean_g1_pathogen_e]] - .data[[mean_g2_pathogen_e]],
        .data[[mean_g1_sterile_e]] - .data[[mean_g2_sterile_e]]
      )
    )
} else {
  df_comparisons <- df_comparisons %>%
    dplyr::mutate(
      diff_tregs_on_minus_off = ifelse(
        injury_type == 'pathogenic',
        .data[[mean_g2_pathogen_e]] - .data[[mean_g1_pathogen_e]],
        .data[[mean_g2_sterile_e]] - .data[[mean_g1_sterile_e]]
      )
    )
}

# Additional difference for randomization check
if (check_randomization_logic) {
  mean_notrnd_pathogen_e <- "mean_tregs_on_pathogen_e"
  mean_notrnd_sterile_e  <- "mean_tregs_on_sterile_e"
  mean_rnd_pathogen_e    <- "mean_tregs_rnd_pathogen_e"
  mean_rnd_sterile_e     <- "mean_tregs_rnd_sterile_e"

  df_comparisons <- df_comparisons %>%
    dplyr::mutate(
      diff_notrnd_minus_rnd = ifelse(
        injury_type == 'pathogenic',
        .data[[mean_notrnd_pathogen_e]] - .data[[mean_rnd_pathogen_e]],
        .data[[mean_notrnd_sterile_e]] - .data[[mean_rnd_sterile_e]]
      )
    )
}

# Select relevant columns
if (check_randomization_logic) {
  df_comparisons <- df_comparisons %>%
    dplyr::select(param_set_id, injury_type, diff_tregs_on_minus_off, diff_notrnd_minus_rnd)
} else {
  df_comparisons <- df_comparisons %>%
    dplyr::select(param_set_id, injury_type, diff_tregs_on_minus_off)
}

# =============================================================================
# MERGE WITH COHEN'S D VALUES
# =============================================================================

# Construct Cohen's d column names
# Try multiple naming conventions
cohens_col_p_v1 <- paste0("d_", group1, "_vs_", group2, "_pathogen_e")
cohens_col_s_v1 <- paste0("d_", group1, "_vs_", group2, "_sterile_e")
cohens_col_p_v2 <- paste0("d_", group2, "_vs_", group1, "_pathogen_e")
cohens_col_s_v2 <- paste0("d_", group2, "_vs_", group1, "_sterile_e")

# Check which naming convention exists
if (cohens_col_p_v1 %in% colnames(df_results)) {
  cohens_col_p <- cohens_col_p_v1
  cohens_col_s <- cohens_col_s_v1
} else if (cohens_col_p_v2 %in% colnames(df_results)) {
  cohens_col_p <- cohens_col_p_v2
  cohens_col_s <- cohens_col_s_v2
} else {
  stop("Error: Could not find Cohen's d columns for ", group1, " vs ", group2, "\n",
       "Tried: ", cohens_col_p_v1, " and ", cohens_col_p_v2)
}

cat("Using Cohen's d columns:\n")
cat("  Pathogenic:", cohens_col_p, "\n")
cat("  Sterile:", cohens_col_s, "\n")

# Merge with Cohen's d
merge_cols <- c('param_set_id', cohens_col_p, cohens_col_s)

if (check_randomization_logic) {
  merge_cols <- c(merge_cols,
                  'd_tregs_rnd_vs_nrnd_pathogen_e',
                  'd_tregs_rnd_vs_nrnd_sterile_e')
}

df_comparisons <- merge(
  df_comparisons,
  distinct(df_results[merge_cols]),
  by = 'param_set_id'
)

# =============================================================================
# CHECK RANDOMIZATION LOGIC (IF APPLICABLE)
# =============================================================================

if (check_randomization_logic) {
  cat("\nChecking randomization logic...\n")

  # Filter based on Cohen's d for tolerance estimation
  df_comparisons_below <- df_comparisons %>%
    dplyr::filter(
      (injury_type == 'pathogenic' & abs(.data[[cohens_col_p]]) < jsd_th) |
      (injury_type == 'sterile' & abs(.data[[cohens_col_s]]) < jsd_th)
    )

  # Dynamic tolerance if not set
  if (is.na(tol_in) || !exists("tol_in")) {
    tol_in <- max(
      max(df_comparisons_below$diff_tregs_on_minus_off),
      abs(min(df_comparisons_below$diff_tregs_on_minus_off))
    )
    cat("  Calculated tol_in:", tol_in, "\n")
  }

  # Check for cases where randomization doesn't make sense
  df_comparisons_keep <- df_comparisons
  df_comparisons_sense <- df_comparisons %>%
    dplyr::mutate(
      makes_sense = ifelse(
        sign(diff_tregs_on_minus_off) != sign(diff_notrnd_minus_rnd) &
        abs(diff_notrnd_minus_rnd) > tol_in &
        (abs(d_tregs_rnd_vs_nrnd_pathogen_e) >= jsd_th |
         abs(d_tregs_rnd_vs_nrnd_sterile_e) >= jsd_th),
        0, 1
      )
    ) %>%
    dplyr::filter(makes_sense == 0)

  cat("  Cases that don't make sense:", nrow(df_comparisons_sense), "\n")

  df_comparisons <- df_comparisons_keep
}

# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

df_comparisons_plot <- df_comparisons

# --- PATHOGENIC ---
df_comparisons_plot_pathogenic <- df_comparisons_plot %>%
  dplyr::filter(injury_type == 'pathogenic')

df_comparisons_plot_pathogenic <- df_comparisons_plot_pathogenic %>%
  dplyr::mutate(
    tregs_better = ifelse(
      diff_tregs_on_minus_off > tol_in, 1,
      ifelse(diff_tregs_on_minus_off < -1 * tol_in, -1, 0)
    )
  )

df_comparisons_plot_pathogenic <- df_comparisons_plot_pathogenic %>%
  dplyr::mutate(
    tregs_better_cohens = ifelse(
      abs(.data[[cohens_col_p]]) > jsd_th,
      tregs_better, 0
    )
  )

# Remove sterile Cohen's d column and rename pathogenic
df_comparisons_plot_pathogenic <- df_comparisons_plot_pathogenic %>%
  dplyr::select(-one_of(cohens_col_s))

colnames(df_comparisons_plot_pathogenic)[
  which(colnames(df_comparisons_plot_pathogenic) == cohens_col_p)
] <- 'cohens_d'

# --- STERILE ---
df_comparisons_plot_sterile <- df_comparisons_plot %>%
  dplyr::filter(injury_type == 'sterile')

df_comparisons_plot_sterile <- df_comparisons_plot_sterile %>%
  dplyr::mutate(
    tregs_better = ifelse(
      diff_tregs_on_minus_off > tol_in, 1,
      ifelse(diff_tregs_on_minus_off < -1 * tol_in, -1, 0)
    )
  )

df_comparisons_plot_sterile <- df_comparisons_plot_sterile %>%
  dplyr::mutate(
    tregs_better_cohens = ifelse(
      abs(.data[[cohens_col_s]]) > jsd_th,
      tregs_better, 0
    )
  )

# Remove pathogenic Cohen's d column and rename sterile
df_comparisons_plot_sterile <- df_comparisons_plot_sterile %>%
  dplyr::select(-one_of(cohens_col_p))

colnames(df_comparisons_plot_sterile)[
  which(colnames(df_comparisons_plot_sterile) == cohens_col_s)
] <- 'cohens_d'

# =============================================================================
# IDENTIFY CONFLICTS BETWEEN STERILE AND PATHOGENIC
# =============================================================================

tregs_better_sterile    <- df_comparisons_plot_sterile %>%
  dplyr::filter(tregs_better_cohens == 1) %>%
  dplyr::pull(param_set_id)

tregs_worse_sterile     <- df_comparisons_plot_sterile %>%
  dplyr::filter(tregs_better_cohens == -1) %>%
  dplyr::pull(param_set_id)

tregs_better_pathogenic <- df_comparisons_plot_pathogenic %>%
  dplyr::filter(tregs_better_cohens == 1) %>%
  dplyr::pull(param_set_id)

tregs_worse_pathogenic  <- df_comparisons_plot_pathogenic %>%
  dplyr::filter(tregs_better_cohens == -1) %>%
  dplyr::pull(param_set_id)

s_better_p_better <- intersect(tregs_better_sterile, tregs_better_pathogenic)
s_better_p_worse  <- intersect(tregs_better_sterile, tregs_worse_pathogenic)
s_worse_p_better  <- intersect(tregs_worse_sterile, tregs_better_pathogenic)
s_worse_p_worse   <- intersect(tregs_worse_sterile, tregs_worse_pathogenic)

cat("\nConflict analysis:\n")
cat("  Better in both:", length(s_better_p_better), "\n")
cat("  Better sterile, worse pathogenic:", length(s_better_p_worse), "\n")
cat("  Worse sterile, better pathogenic:", length(s_worse_p_better), "\n")
cat("  Worse in both:", length(s_worse_p_worse), "\n")

# =============================================================================
# COMBINE AND SAVE RESULTS
# =============================================================================

df_comparisons_plot <- rbind(df_comparisons_plot_pathogenic, df_comparisons_plot_sterile)
output_rds <- paste0('./df_comparisons_plot_', output_prefix, '.rds')
saveRDS(df_comparisons_plot, output_rds)
cat("\nSaved results to:", output_rds, "\n")

# =============================================================================
# GENERATE PLOT
# =============================================================================

cohens_th <- jsd_th
e_th      <- tol_in

source('./MISC/REGIONS.R')

output_png <- paste0("./ABM_JSD_", toupper(output_prefix), ".png")

ggsave(
  filename = output_png,
  plot = p,
  width = 9,
  height = 6,
  dpi = 300,
  bg = 'white'
)

cat("Saved plot to:", output_png, "\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total unique param sets:", length(unique(df_comparisons_plot$param_set_id)), "\n")

# Merge pathogenic and sterile for summary
dfp_p <- df_comparisons_plot %>%
  dplyr::filter(injury_type == 'pathogenic') %>%
  dplyr::select(param_set_id, tregs_better_cohens)

dfp_s <- df_comparisons_plot %>%
  dplyr::filter(injury_type == 'sterile') %>%
  dplyr::select(param_set_id, tregs_better_cohens)

colnames(dfp_p)[2] <- 'tregs_better_cohens_p'
colnames(dfp_s)[2] <- 'tregs_better_cohens_s'

dfp_sp <- merge(dfp_p, dfp_s, by = 'param_set_id')
dfp_sp <- dfp_sp %>%
  dplyr::mutate(
    treg_better_cohens_both = ifelse(tregs_better_cohens_p + tregs_better_cohens_s == 2, 1, 0)
  )

cat("\nPathogenic injury:\n")
cat("  Better:",
    round(100 * length(which(dfp_sp$tregs_better_cohens_p == 1)) / nrow(dfp_sp), 1),
    "%\n")
cat("  Worse:",
    round(100 * length(which(dfp_sp$tregs_better_cohens_p == -1)) / nrow(dfp_sp), 1),
    "%\n")

cat("\nSterile injury:\n")
cat("  Better:",
    round(100 * length(which(dfp_sp$tregs_better_cohens_s == 1)) / nrow(dfp_sp), 1),
    "%\n")
cat("  Worse:",
    round(100 * length(which(dfp_sp$tregs_better_cohens_s == -1)) / nrow(dfp_sp), 1),
    "%\n")

cat("\nBetter in both injuries:",
    round(100 * sum(dfp_sp$treg_better_cohens_both) / nrow(dfp_sp), 1),
    "%\n")

# Non-zero cases
df_comparisons_plot_nz <- df_comparisons_plot %>%
  dplyr::filter(tregs_better_cohens != 0)

cat("\nNon-zero cases:", nrow(df_comparisons_plot_nz), "\n")

# Print conflicts
cat("\nConflicting param sets (better in one, worse in other):\n")
conflicts <- c(s_better_p_worse, s_worse_p_better)
if (length(conflicts) > 0) {
  print(conflicts)
} else {
  cat("  None\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
