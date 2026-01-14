# ============================================================================
# TEST SCRIPT FOR NEW C++ OPTIMIZATIONS
# ============================================================================
# This script tests the compilation and basic functionality of the new C++
# functions in FAST_FUNCTIONS.cpp
#
# Run this script BEFORE using MAIN_SIM_FASTER.R to ensure all C++ functions
# compile correctly.
# ============================================================================

cat("Testing FAST_FUNCTIONS.cpp compilation...\n\n")

# Compile C++ functions
cat("Step 1: Compiling C++ code...\n")
Rcpp::sourceCpp('MISC/FAST_FUNCTIONS.cpp')
cat("✓ Compilation successful!\n\n")

# Test basic functionality
cat("Step 2: Testing new C++ functions...\n")

# Test 1: move_cells_chemotaxis_cpp
cat("  Testing move_cells_chemotaxis_cpp... ")
test_x = c(5, 10, 15)
test_y = c(5, 10, 15)
test_density = matrix(runif(400), 20, 20)
result = move_cells_chemotaxis_cpp(test_x, test_y, test_density, 20)
stopifnot(length(result$x) == 3)
stopifnot(length(result$y) == 3)
cat("✓\n")

# Test 2: update_PAMPs_cpp
cat("  Testing update_PAMPs_cpp... ")
test_PAMPs = matrix(0, 20, 20)
test_pathogen_coords = matrix(c(5, 5, 1, 10, 10, 2), ncol=3, byrow=TRUE)
colnames(test_pathogen_coords) = c("x", "y", "id")
result = update_PAMPs_cpp(test_PAMPs, test_pathogen_coords, 1.0, 1.0, 5.0, 5.0)
stopifnot(nrow(result) == 20)
stopifnot(ncol(result) == 20)
cat("✓\n")

# Test 3: process_engulfment_batch_cpp
cat("  Testing process_engulfment_batch_cpp... ")
test_phago_x = c(5, 10)
test_phago_y = c(5, 10)
test_activity = c(0.5, 0.8)
test_phenotype = c(0L, 1L)
test_pathogen_coords = matrix(c(5, 5, 1, 10, 10, 2), ncol=3, byrow=TRUE)
colnames(test_pathogen_coords) = c("x", "y", "id")
test_commensal_coords = matrix(c(5, 5, 3), ncol=3, byrow=TRUE)
colnames(test_commensal_coords) = c("x", "y", "id")
test_registry = matrix(0, nrow=2, ncol=10)
result = process_engulfment_batch_cpp(
  test_phago_x, test_phago_y, test_activity, test_phenotype,
  test_pathogen_coords, test_commensal_coords, test_registry, 10L
)
stopifnot(!is.null(result$pathogen_coords))
stopifnot(!is.null(result$commensal_coords))
cat("✓\n")

# Test 4: activate_tregs_batch_cpp
cat("  Testing activate_tregs_batch_cpp... ")
test_M_indices = c(1L, 2L)
test_phago_x = c(5, 10)
test_phago_y = c(5, 10)
test_pat_engulfed = c(2L, 3L)
test_com_engulfed = c(1L, 2L)
test_treg_x = c(6, 11)
test_treg_y = c(6, 11)
test_treg_pheno = c(0L, 0L)
test_treg_activity = c(0, 0)
test_treg_age = c(0L, 0L)
result = activate_tregs_batch_cpp(
  test_M_indices, test_phago_x, test_phago_y,
  test_pat_engulfed, test_com_engulfed,
  test_treg_x, test_treg_y, test_treg_pheno,
  test_treg_activity, test_treg_age, 2L, 0.9
)
stopifnot(length(result$treg_phenotype) == 2)
cat("✓\n")

# Test 5: update_M0_phenotypes_cpp
cat("  Testing update_M0_phenotypes_cpp... ")
test_M0_indices = c(1L, 2L)
test_phago_x = c(5, 10)
test_phago_y = c(5, 10)
test_DAMPs = matrix(runif(400), 20, 20)
test_SAMPs = matrix(runif(400), 20, 20)
test_PAMPs = matrix(runif(400), 20, 20)
result = update_M0_phenotypes_cpp(
  test_M0_indices, test_phago_x, test_phago_y,
  test_DAMPs, test_SAMPs, test_PAMPs,
  1L, 1L, 1L, 20L, 0.5, 0.5, 0.1, 0.5
)
stopifnot(length(result$phenotype) == 2)
cat("✓\n")

# Test 6: update_active_phenotypes_cpp
cat("  Testing update_active_phenotypes_cpp... ")
test_active_indices = c(1L, 2L)
test_phago_x = c(5, 10)
test_phago_y = c(5, 10)
test_pheno = c(1L, 2L)
test_age = c(5L, 5L)
result = update_active_phenotypes_cpp(
  test_active_indices, test_phago_x, test_phago_y,
  test_pheno, test_age,
  test_DAMPs, test_SAMPs, test_PAMPs,
  1L, 1L, 1L, 20L, 10L, 0.5, 0.5,
  0.05, 0.3, 0.1, 0.5, 0.05, 0.4,
  only_suppress = FALSE
)
stopifnot(length(result$phenotype) == 2)
cat("✓\n")

# Test 7: shift_registry_batch_cpp
cat("  Testing shift_registry_batch_cpp... ")
test_registry = matrix(c(1, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE)
test_should_shift = c(TRUE, FALSE)
result = shift_registry_batch_cpp(test_registry, test_should_shift)
stopifnot(nrow(result) == 2)
stopifnot(ncol(result) == 3)
stopifnot(result[1, 1] == 0)  # First row should be shifted
stopifnot(result[2, 1] == 4)  # Second row should be unchanged
cat("✓\n")

cat("\n")
cat("========================================\n")
cat("ALL TESTS PASSED! ✓\n")
cat("========================================\n")
cat("\n")
cat("The new C++ functions are working correctly.\n")
cat("You can now use MAIN_SIM_FASTER.R for improved performance.\n")
