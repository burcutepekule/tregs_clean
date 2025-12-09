# Configuration file for running different comparisons
# This script sets up parameters and sources the generic comparison script

# ============= COMMON PARAMETERS =============
# These can be overridden for specific comparisons if needed
jsd_th = 0.3
tol_in = 125 * 0.2  # Can be set to NA for auto-calculation
filter_control = 0

# ============= FUNCTION TO RUN A COMPARISON =============
run_comparison = function(comparison_name,
                          group1,
                          group2,
                          jsd_th = 0.3,
                          tol_in = NA,
                          filter_control = 0,
                          perform_randomness_check = FALSE,
                          perform_macspec_analysis = FALSE) {

  # Assign to global environment so the sourced script can access them
  assign("comparison_name", comparison_name, envir = .GlobalEnv)
  assign("group1", group1, envir = .GlobalEnv)
  assign("group2", group2, envir = .GlobalEnv)
  assign("jsd_th", jsd_th, envir = .GlobalEnv)
  assign("tol_in", tol_in, envir = .GlobalEnv)
  assign("filter_control", filter_control, envir = .GlobalEnv)
  assign("perform_randomness_check", perform_randomness_check, envir = .GlobalEnv)
  assign("perform_macspec_analysis", perform_macspec_analysis, envir = .GlobalEnv)

  cat("\n", strrep("=", 80), "\n")
  cat("Running comparison:", comparison_name, "\n")
  cat("Group 1:", group1, "vs Group 2:", group2, "\n")
  cat(strrep("=", 80), "\n\n")

  # Source the generic script
  source('./DLL_datazanalyse_abm_regions_generic.R')

  cat("\n", strrep("=", 80), "\n")
  cat("Completed comparison:", comparison_name, "\n")
  cat(strrep("=", 80), "\n\n")
}

# ============= PREDEFINED COMPARISONS =============

# Comparison 1: Tregs ON vs Tregs OFF
run_comparison_tregs_on_vs_off = function() {
  run_comparison(
    comparison_name = "tregs_on_vs_tregs_off",
    group1 = "tregs_on",
    group2 = "tregs_off",
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0,
    perform_randomness_check = TRUE,
    perform_macspec_analysis = FALSE
  )
}

# Comparison 2: Tregs ON vs Macspec
run_comparison_tregs_on_vs_macspec = function() {
  run_comparison(
    comparison_name = "tregs_on_vs_macspec",
    group1 = "tregs_on",
    group2 = "macspec",
    jsd_th = 0.3,
    tol_in = NA,  # Auto-calculate
    filter_control = 1,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = TRUE
  )
}

# Comparison 3: Macspec vs Tregs OFF
run_comparison_macspec_vs_tregs_off = function() {
  run_comparison(
    comparison_name = "macspec_vs_tregs_off",
    group1 = "macspec",
    group2 = "tregs_off",
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = FALSE
  )
}

# ============= ADDITIONAL COMPARISON EXAMPLES =============

# Comparison 4: Control vs Tregs OFF
run_comparison_ctrl_vs_tregs_off = function() {
  run_comparison(
    comparison_name = "ctrl_vs_tregs_off",
    group1 = "ctrl",
    group2 = "tregs_off",
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = FALSE
  )
}

# Comparison 5: Tregs RND vs Tregs OFF
run_comparison_tregs_rnd_vs_tregs_off = function() {
  run_comparison(
    comparison_name = "tregs_rnd_vs_tregs_off",
    group1 = "tregs_rnd",
    group2 = "tregs_off",
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = FALSE
  )
}

# Comparison 6: Tregs ON vs Tregs RND
run_comparison_tregs_on_vs_tregs_rnd = function() {
  run_comparison(
    comparison_name = "tregs_on_vs_tregs_rnd",
    group1 = "tregs_on",
    group2 = "tregs_rnd",
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = FALSE
  )
}

# ============= RUN ALL COMPARISONS =============
run_all_comparisons = function() {
  run_comparison_tregs_on_vs_off()
  run_comparison_tregs_on_vs_macspec()
  run_comparison_macspec_vs_tregs_off()
}

# ============= INSTRUCTIONS =============
cat("
===================================================================
Generic Comparison Script Configuration
===================================================================

Available functions:
  - run_comparison_tregs_on_vs_off()
  - run_comparison_tregs_on_vs_macspec()
  - run_comparison_macspec_vs_tregs_off()
  - run_comparison_ctrl_vs_tregs_off()
  - run_comparison_tregs_rnd_vs_tregs_off()
  - run_comparison_tregs_on_vs_tregs_rnd()
  - run_all_comparisons()

Custom comparisons:
  run_comparison(
    comparison_name = 'my_comparison',
    group1 = 'group1_name',
    group2 = 'group2_name',
    jsd_th = 0.3,
    tol_in = NA,
    filter_control = 0,
    perform_randomness_check = FALSE,
    perform_macspec_analysis = FALSE
  )

Examples:
  # Run a single comparison
  run_comparison_tregs_on_vs_off()

  # Run all predefined comparisons
  run_all_comparisons()

  # Custom comparison
  run_comparison('custom', 'tregs_on', 'ctrl', jsd_th=0.5)

===================================================================
")
