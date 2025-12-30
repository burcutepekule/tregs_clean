rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)
library(ggsignif)

jsd_th         = 0.3
tol_in_e       = 125*0.25
tol_in_p       = -1*25*25*0.05
# tol_in_e       = 0
# tol_in_p       = 0
M1_M2_diff     = 0
filter_control = 0
labels_on      = 1
score_type     = 'pathogen'
data_suffix    = '_ros_vs_ctrl_nopatlevel' #
data_suffix    = '_ros_vs_ctrl_patros' #
data_suffix    = '_patros' #

source("./MISC/PLOT_FUNCTIONS_ABM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('./lhs_parameters_della.csv', show_col_types = FALSE)
df_results_keep = readRDS(paste0('./data_cpp_read_abm',data_suffix,'.rds'))
unique(df_results_keep$param_set_id)
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df = as.data.frame(table(df_results_keep$param_set_id))
reps_df$Var1 = as.numeric(as.character(reps_df$Var1))
if(data_suffix == '_ros_vs_ctrl_nopatlevel'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==40) %>% dplyr::pull(Var1) # 40 = 10 reps per scenario, 2 scenarios x 2 times recording for epithelial and pathogen scores 
}else if(data_suffix == '_ros_vs_ctrl_patros'){
  keep_param_id = reps_df %>% dplyr::filter(Freq==260) %>% dplyr::pull(Var1) # 260 = 10 reps per scenario, 13 scenarios x 2 times recording for epithelial and pathogen scores 
}else if(data_suffix == '_patros'){
  # keep_param_id = reps_df %>% dplyr::filter(Freq==1100) %>% dplyr::pull(Var1) # 1100 = 5 reps per scenario, 110 scenarios x 2 times recording for epithelial and pathogen scores
  keep_param_id = reps_df %>% dplyr::filter(Freq > 220) %>% dplyr::pull(Var1) # 1100 = 5 reps per scenario, 110 scenarios x 2 times recording for epithelial and pathogen scores
  # keep_param_id = reps_df %>% dplyr::pull(Var1)
}

df_results = df_results_keep %>% filter(param_set_id %in% keep_param_id)
length(unique(df_results$param_set_id))

### ----- filter based on ss_start, it cannot be too large otherwise not much to compare!
# ss_start_threshold = 4500 # used to be 4500, just for simulation purposes to save time
ss_start_threshold = 1800 
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
num_params = length(unique(df_results$param_set_id))
print(num_params)
saveRDS(unique(df_results$param_set_id), 'evo_selected_level_0.rds')

df_comparisons = df_results %>% dplyr::select(
  param_set_id, injury_type, 
  starts_with("d_ros"), starts_with("mean_ros"))

df_comparisons_keep = distinct(df_comparisons)

# df_comparisons_keep_pick = df_comparisons_keep %>% dplyr::filter(param_set_id==4201) %>% dplyr::select(-injury_type)
# 
# df_comparisons_keep_pick = pivot_longer(
#   df_comparisons_keep_pick, 
#   cols = -c(param_set_id),  # or specify which cols TO pivot
#   names_to = 'variable',
#   values_to = 'value'
# )

# ============= HELPER FUNCTION: Define "infection under control" ====================================================
# Infection is under control when:
# 1. Epithelial health is approximately 150 (threshold: > 149)
# 2. Pathogen load is approximately 0 (threshold: < 1)
is_under_control = function(epithelial_health, pathogen_load) {
  return(epithelial_health > 150*0.75 & pathogen_load < 10)
}

# ============= STEP 1: Filter for evolutionary need for ROS ====================================================
# Rule 1: Without ROS (ros_0, treg0_), even the lowest pathogen level (pat1_) should:
#   a) Cause significant injury to epithelium (ros0_pat1 NOT under control)
#   b) Have higher pathogen burden than ros1_pat1 case (ros1_pat1 should be better)

# Start with the base filtering
df_steps = df_comparisons_keep %>%
  dplyr::filter(injury_type == 'pathogenic')

ros_vals = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
pat_vals = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Dynamically create the controlled status columns for all ros/pat combinations
for (ros in ros_vals) {
  for (pat in pat_vals) {
    col_name = paste0("ros", ros, "_pat", pat, "_controlled")
    mean_e_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_e")
    mean_p_col = paste0("mean_ros", ros, "_pat", pat, "_treg0_p")
    
    df_steps = df_steps %>%
      dplyr::mutate(!!col_name := is_under_control(.data[[mean_e_col]], .data[[mean_p_col]]))
  }
}

# # ==== For each parameter set, check if it has the triangular pattern

# Parameter	Default	What it controls	Adjust to...
# min_first_row_ratio (0.7)  : How wide the funnel top is - Increase (0.8-0.9) for wider top, decrease (0.5-0.6) for narrower top
# n_bottom_rows (3)          : How many bottom rows to check 2-5 rows typical
# max_bottom_ratio (0.15)    : How closed the funnel bottom is - Decrease (0.05-0.10) for tighter closure, increase (0.20-0.30) for partial closure
# min_top_bottom_ratio (1.5) : Strength of funnel taper - Increase (2.0-3.0) for dramatic narrowing, decrease (1.2-1.3) for gentle taper
# max_center_row_ratio (0.4) : How top-heavy the pattern is - Decrease (0.3) for more top-heavy, increase (0.5-0.6) for more distributed
# min_total_ratio	(0.15)     : Minimum density of TRUEs - Increase to require more controlled regions
# max_total_ratio	(0.65)     : Maximum density of TRUEs - Decrease to require more selectivity

# # ==== To find MORE matches (more permissive):
#   
# # min_first_row_ratio = 0.6      # Allow narrower tops
# # max_bottom_ratio = 0.25        # Allow more open bottoms
# # min_top_bottom_ratio = 1.2     # Weaker taper OK
# # max_total_ratio = 0.75         # Allow more TRUEs overall
# 
# # ==== To find FEWER, stricter matches:
# 
# # min_first_row_ratio = 0.8      # Require wider tops
# # max_bottom_ratio = 0.10        # Require tighter closure
# # min_top_bottom_ratio = 2.0     # Require stronger taper
# # max_total_ratio = 0.55         # Limit total TRUEs

# Detect irregular triangular/blob pattern
# Detect inverted funnel pattern (wide at top, narrow/closed at bottom)
df_steps$has_triangular_pattern = apply(df_steps, 1, function(row) {
  
  # Build control matrix
  control_matrix = matrix(FALSE, nrow = length(pat_vals), ncol = length(ros_vals))
  for (i in seq_along(pat_vals)) {
    for (j in seq_along(ros_vals)) {
      pat = pat_vals[i]
      ros = ros_vals[j]
      col_name = paste0("ros", ros, "_pat", pat, "_controlled")
      control_matrix[i, j] = as.logical(row[[col_name]])
    }
  }
  
  n_true_per_row = rowSums(control_matrix)
  
  # ========== TUNABLE PARAMETERS ==========
  
  # PARAMETER 1: How "open" should the top be?
  min_first_row_ratio = 0.7  # 0.7 = 70% of first row must be TRUE
  # Increase (e.g., 0.8, 0.9) for wider top
  # Decrease (e.g., 0.5, 0.6) for narrower top
  
  first_row_ratio = n_true_per_row[1] / ncol(control_matrix)
  if (first_row_ratio < min_first_row_ratio) return(FALSE)
  
  
  # PARAMETER 2: How "closed" should the bottom be?
  n_bottom_rows = 3          # How many rows from bottom to check (3-5 typical)
  max_bottom_ratio = 0.15    # Maximum 15% of bottom rows can be TRUE
  # Decrease (e.g., 0.05, 0.10) for tighter closure
  # Increase (e.g., 0.20, 0.30) for looser closure
  
  last_n_rows = tail(n_true_per_row, n_bottom_rows)
  last_rows_ratio = sum(last_n_rows) / (n_bottom_rows * ncol(control_matrix))
  if (last_rows_ratio > max_bottom_ratio) return(FALSE)
  
  
  # PARAMETER 3: How strong should the narrowing trend be?
  # (Compares top half vs bottom half)
  min_top_bottom_ratio = 1.5  # Top half should have 1.5x more TRUEs than bottom
  # Increase (e.g., 2.0, 2.5) for stronger funnel
  # Decrease (e.g., 1.2, 1.3) for gentler taper
  
  mid_point = ceiling(length(n_true_per_row) / 2)
  top_half_mean = mean(n_true_per_row[1:mid_point])
  bottom_half_mean = mean(n_true_per_row[(mid_point+1):length(n_true_per_row)])
  
  if (bottom_half_mean == 0) {
    # Bottom is completely empty - that's good for funnel
    top_bottom_ratio = Inf
  } else {
    top_bottom_ratio = top_half_mean / bottom_half_mean
  }
  
  if (top_bottom_ratio < min_top_bottom_ratio) return(FALSE)
  
  
  # PARAMETER 4: Where should the "center of mass" be?
  max_center_row_ratio = 0.4  # Center should be in top 40% of rows
  # Decrease (e.g., 0.3) for more top-heavy
  # Increase (e.g., 0.5) for more centered
  
  if (sum(control_matrix) >= 3) {  # Need at least some TRUEs
    true_positions = which(control_matrix, arr.ind = TRUE)
    center_row = mean(true_positions[, 1])
    
    if (center_row > length(pat_vals) * max_center_row_ratio) return(FALSE)
  } else {
    return(FALSE)  # Too few TRUEs
  }
  
  
  # PARAMETER 5: How much total coverage?
  min_total_ratio = 0.15  # At least 15% of all cells should be TRUE
  max_total_ratio = 0.65  # At most 65% of all cells should be TRUE
  # Narrow this range for stricter funnel shape
  # Widen for more permissive matching
  
  total_true_ratio = sum(control_matrix) / length(control_matrix)
  if (total_true_ratio < min_total_ratio || total_true_ratio > max_total_ratio) {
    return(FALSE)
  }
  
  return(TRUE)
})

# Find param_set_ids with triangular pattern
triangular_pattern_ids = df_steps$param_set_id[df_steps$has_triangular_pattern == TRUE]

cat("Found", length(triangular_pattern_ids), "parameter sets with triangular pattern\n")
print(triangular_pattern_ids)

# Save them
saveRDS(triangular_pattern_ids, 'evo_selected_triangular_pattern.rds')

# Find param_set_ids with the CORRECT triangular pattern
triangular_pattern_ids = df_steps$param_set_id[df_steps$has_triangular_pattern == TRUE]

cat("Found", length(triangular_pattern_ids), "parameter sets with correct triangular pattern\n")
print(triangular_pattern_ids)

# Visualize to verify
visualize_control_pattern = function(param_id) {
  row = df_steps[df_steps$param_set_id == param_id, ]
  
  control_matrix = matrix(NA, nrow = length(pat_vals), ncol = length(ros_vals))
  rownames(control_matrix) = paste0("pat", pat_vals)
  colnames(control_matrix) = paste0("ros", ros_vals)
  
  for (i in seq_along(pat_vals)) {
    for (j in seq_along(ros_vals)) {
      pat = pat_vals[i]
      ros = ros_vals[j]
      col_name = paste0("ros", ros, "_pat", pat, "_controlled")
      control_matrix[i, j] = as.logical(row[[col_name]])
    }
  }
  
  cat("\nControl pattern for param_set_id", param_id, "\n")
  cat("(TRUE = controlled, FALSE = not controlled)\n\n")
  print(control_matrix)
  
  # Calculate min ROS needed for each pat level
  min_ros_needed = sapply(1:nrow(control_matrix), function(i) {
    true_indices = which(control_matrix[i, ])
    if (length(true_indices) > 0) ros_vals[min(true_indices)] else NA
  })
  
  cat("\nMinimum ROS needed for each pathogen level:\n")
  result_df = data.frame(pat_level = pat_vals, min_ros_needed = min_ros_needed)
  print(result_df)
  
  # Show if pattern is increasing
  cat("\nPattern check: min ROS should increase with pat level\n")
  cat("Differences between consecutive rows:", diff(min_ros_needed[!is.na(min_ros_needed)]), "\n")
}

# Test with first result
if (length(triangular_pattern_ids) > 0) {
  visualize_control_pattern(triangular_pattern_ids[1])
}

# # Save all control matrices in a single list
# control_matrices_list = list()
# 
# for (param_id in triangular_pattern_ids) {
#   row = df_steps[df_steps$param_set_id == param_id, ]
#   
#   # Build the control matrix
#   control_matrix = matrix(NA, nrow = length(pat_vals), ncol = length(ros_vals))
#   rownames(control_matrix) = paste0("pat", pat_vals)
#   colnames(control_matrix) = paste0("ros", ros_vals)
#   
#   for (i in seq_along(pat_vals)) {
#     for (j in seq_along(ros_vals)) {
#       pat = pat_vals[i]
#       ros = ros_vals[j]
#       col_name = paste0("ros", ros, "_pat", pat, "_controlled")
#       control_matrix[i, j] = as.logical(row[[col_name]])
#     }
#   }
#   
#   # Add to list with param_id as name
#   control_matrices_list[[as.character(param_id)]] = control_matrix
# }
# 
# # Save the entire list
# saveRDS(control_matrices_list, 'control_matrices_all_triangular_patterns.rds')
# 
# # To load and access later:
# # matrices = readRDS('control_matrices_all_triangular_patterns.rds')
# # matrices[["20601"]]  # Access specific param_set_id
# 
# # Save
# saveRDS(triangular_pattern_ids, 'evo_selected_triangular_pattern.rds')
