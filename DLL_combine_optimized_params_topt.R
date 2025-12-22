rm(list=ls())
library(dplyr)
library(tidyr)

# ============================================================================
# COMBINE OPTIMIZED TREG PARAMETERS WITH BASE PARAMETERS
# ============================================================================

cat("Loading base parameters...\n")
params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
treg_param_bounds = readRDS("./treg_param_bounds.rds")

cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# LOAD ALL OPTIMIZATION RESULTS
# ============================================================================

optimization_dir = "./optimization_results"
if (!dir.exists(optimization_dir)) {
  stop("Optimization results directory not found: ", optimization_dir)
}

opt_files = list.files(optimization_dir, pattern = "optimized_treg_params_.*\\.rds", full.names = TRUE)
cat("Found", length(opt_files), "optimization result files\n\n")

if (length(opt_files) == 0) {
  stop("No optimization results found!")
}

# Load all optimization results
all_opt_results = list()
for (file in opt_files) {
  result = readRDS(file)
  key = paste0("set_", result$param_set_id, "_scenario_", result$scenario_ind)
  all_opt_results[[key]] = result
}

cat("Loaded", length(all_opt_results), "optimization results\n\n")

# ============================================================================
# CREATE FULL PARAMETER SETS WITH OPTIMIZED TREG PARAMETERS
# ============================================================================

# For each unique param_set_id and scenario combination, create a full parameter set
full_params_list = list()

for (key in names(all_opt_results)) {
  result = all_opt_results[[key]]

  param_set_id = result$param_set_id
  scenario_ind = result$scenario_ind

  # Get base parameters
  base_params = params_df %>% filter(param_set_id == !!param_set_id)

  if (nrow(base_params) == 0) {
    cat("Warning: No base parameters found for param_set_id", param_set_id, "\n")
    next
  }

  # Add optimized Treg parameters
  full_params = base_params
  for (param_name in names(result$optimized_params)) {
    full_params[[param_name]] = result$optimized_params[param_name]
  }

  # Add scenario information
  full_params$scenario_ind = scenario_ind
  full_params$optimization_converged = result$convergence == 0
  full_params$objective_value = result$objective_value

  full_params_list[[key]] = full_params
}

# Combine all parameter sets
full_params_df = bind_rows(full_params_list)

cat("Created", nrow(full_params_df), "full parameter sets with optimized Treg parameters\n\n")

# ============================================================================
# EXPORT COMBINED PARAMETERS
# ============================================================================

# Save full parameter sets
output_file = "./lhs_parameters_with_optimized_tregs.csv"
write.csv(full_params_df, output_file, row.names = FALSE)
cat("Full parameter sets saved to:", output_file, "\n")

# Also create a summary of optimized Treg parameters
summary_df = full_params_df %>%
  select(param_set_id, scenario_ind,
         rat_com_pat_threshold, diffusion_speed_SAMPs, add_SAMPs,
         SAMPs_decay, activation_threshold_SAMPs, treg_discrimination_efficiency,
         optimization_converged, objective_value)

summary_file = "./optimized_treg_params_summary.csv"
write.csv(summary_df, summary_file, row.names = FALSE)
cat("Summary of optimized Treg parameters saved to:", summary_file, "\n")

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Converged optimizations:", sum(full_params_df$optimization_converged), "/", nrow(full_params_df), "\n\n")

cat("Optimized Treg parameter ranges:\n")
for (param_name in names(treg_param_bounds)) {
  if (param_name %in% names(full_params_df)) {
    values = full_params_df[[param_name]]
    cat(sprintf("  %s: [%.3f, %.3f] (mean: %.3f)\n",
                param_name, min(values), max(values), mean(values)))
  }
}

cat("\nDone!\n")
