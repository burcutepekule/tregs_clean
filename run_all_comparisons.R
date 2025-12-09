# =============================================================================
# Run All Comparisons Script
# =============================================================================
# This script runs all three original analyses using the generic script
# =============================================================================

# Set common parameters
jsd_th <- 0.3    # Cohen's d threshold
tol_in <- 18.75  # Tolerance for difference

# =============================================================================
# Analysis 1: macspec vs tregs_off (Baseline comparison)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS 1: macspec vs tregs_off\n")
cat("================================================================================\n")

group1 <- "macspec"
group2 <- "tregs_off"
difference_direction <- "group1_minus_group2"
filter_control <- 0
output_prefix <- "macspec_vs_tregs_off"
check_randomization_logic <- FALSE
show_intermediate_picks <- FALSE

source("./DLL_datazanalyse_abm_regions_generic.R")

# =============================================================================
# Analysis 2: tregs_on vs macspec
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS 2: tregs_on vs macspec\n")
cat("================================================================================\n")

group1 <- "tregs_on"
group2 <- "macspec"
difference_direction <- "group1_minus_group2"
filter_control <- 1
output_prefix <- "tregs_on_vs_macspec"
check_randomization_logic <- FALSE
show_intermediate_picks <- TRUE

source("./DLL_datazanalyse_abm_regions_generic.R")

# =============================================================================
# Analysis 3: tregs_on vs tregs_off
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS 3: tregs_on vs tregs_off\n")
cat("================================================================================\n")

group1 <- "tregs_on"
group2 <- "tregs_off"
difference_direction <- "group1_minus_group2"
filter_control <- 1
output_prefix <- "tregs_on_vs_tregs_off"
check_randomization_logic <- TRUE
show_intermediate_picks <- FALSE

source("./DLL_datazanalyse_abm_regions_generic.R")

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("================================================================================\n")
cat("Generated files:\n")
cat("  - df_comparisons_plot_macspec_vs_tregs_off.rds\n")
cat("  - ABM_JSD_MACSPEC_VS_TREGS_OFF.png\n")
cat("  - df_comparisons_plot_tregs_on_vs_macspec.rds\n")
cat("  - ABM_JSD_TREGS_ON_VS_MACSPEC.png\n")
cat("  - df_comparisons_plot_tregs_on_vs_tregs_off.rds\n")
cat("  - ABM_JSD_TREGS_ON_VS_TREGS_OFF.png\n")
cat("================================================================================\n")
