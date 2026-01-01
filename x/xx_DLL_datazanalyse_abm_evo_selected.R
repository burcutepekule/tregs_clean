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

# For each parameter set, check if it has the triangular pattern

# Checks first row has many TRUEs (≥60%) - low pathogen should be controllable
# Uses correlation to verify decreasing trend in number of TRUEs per row
# Verifies diagonal shift - first TRUE position increases with pat level
# Requires bottom rows to be uncontrollable - high pathogen can't be controlled
# Checks total TRUE percentage - avoids edge cases

df_steps$has_triangular_pattern = apply(df_steps, 1, function(row) {
  
  # Build a matrix: rows = pat_vals, cols = ros_vals
  control_matrix = matrix(NA, nrow = length(pat_vals), ncol = length(ros_vals))
  
  for (i in seq_along(pat_vals)) {
    for (j in seq_along(ros_vals)) {
      pat = pat_vals[i]
      ros = ros_vals[j]
      col_name = paste0("ros", ros, "_pat", pat, "_controlled")
      control_matrix[i, j] = as.logical(row[[col_name]])
    }
  }
  
  # Calculate key metrics for each row
  n_true_per_row = rowSums(control_matrix, na.rm = TRUE)
  
  # For each row, find first and last TRUE position
  first_true_idx = sapply(1:nrow(control_matrix), function(i) {
    true_indices = which(control_matrix[i, ])
    if (length(true_indices) > 0) min(true_indices) else Inf
  })
  
  last_true_idx = sapply(1:nrow(control_matrix), function(i) {
    true_indices = which(control_matrix[i, ])
    if (length(true_indices) > 0) max(true_indices) else -Inf
  })
  
  # ===== TRIANGULAR PATTERN REQUIREMENTS =====
  
  # 1. First row (lowest pat) should be mostly controlled (>50% TRUEs)
  if (n_true_per_row[1] < length(ros_vals) * 0.5) {
    return(FALSE)
  }
  
  # 2. Number of TRUEs should decrease as pat level increases
  # (using correlation - should be negative)
  n_true_trend = cor(1:length(n_true_per_row), n_true_per_row, use = "complete.obs")
  if (is.na(n_true_trend) || n_true_trend > -0.4) {
    return(FALSE)
  }
  
  # 3. For controllable rows, first TRUE should shift right (need more ROS)
  controllable_rows = which(is.finite(first_true_idx))
  if (length(controllable_rows) >= 3) {
    first_true_controllable = first_true_idx[controllable_rows]
    
    # Check if mostly increasing (allowing some noise)
    diffs = diff(first_true_controllable)
    prop_nondecreasing = sum(diffs >= 0) / length(diffs)
    
    if (prop_nondecreasing < 0.6) {
      return(FALSE)
    }
  }
  
  # 4. Last rows (high pat) should be mostly uncontrollable
  # At least 50% of last 4 rows should be all FALSE
  n_bottom_rows = min(4, length(pat_vals))
  bottom_rows = (length(pat_vals) - n_bottom_rows + 1):length(pat_vals)
  n_uncontrollable_bottom = sum(n_true_per_row[bottom_rows] == 0)
  
  if (n_uncontrollable_bottom < n_bottom_rows * 0.5) {
    return(FALSE)
  }
  
  # 5. Total number of TRUEs should be reasonable (not too few, not all)
  total_true = sum(n_true_per_row)
  total_cells = length(pat_vals) * length(ros_vals)
  
  if (total_true < total_cells * 0.15 || total_true > total_cells * 0.8) {
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
