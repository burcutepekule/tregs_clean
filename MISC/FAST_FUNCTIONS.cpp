// ============================================================================
// HIGH-PERFORMANCE C++ IMPLEMENTATIONS FOR TREGS SIMULATION
// ============================================================================
// This file contains C++ implementations of computational bottlenecks
// Expected speedup: 20-100x faster than R for these operations
//
// Compile with: Rcpp::sourceCpp("MISC/FAST_FUNCTIONS.cpp")
// ============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// ============================================================================
// SIGNAL CALCULATIONS (MAJOR BOTTLENECK)
// ============================================================================

// [[Rcpp::export]]
double get_8n_avg_signal_cpp(int x, int y, int act_radius_signal,
                             NumericMatrix signal_matrix, int grid_size)
{
  // Convert from R's 1-based to C++'s 0-based indexing
  int x0 = x - 1;
  int y0 = y - 1;

  // Calculate bounds with edge checking
  int x_start = std::max(0, x0 - act_radius_signal);
  int x_end = std::min(grid_size - 1, x0 + act_radius_signal);
  int y_start = std::max(0, y0 - act_radius_signal);
  int y_end = std::min(grid_size - 1, y0 + act_radius_signal);

  double sum = 0.0;
  int count = 0;

  for (int yi = y_start; yi <= y_end; yi++)
  {
    for (int xi = x_start; xi <= x_end; xi++)
    {
      sum += signal_matrix(yi, xi);
      count++;
    }
  }

  return (count > 0) ? (sum / count) : 0.0;
}

// Vectorized version - process multiple agents at once
// [[Rcpp::export]]
NumericVector get_8n_avg_signal_vectorized_cpp(
    IntegerVector x_vec, IntegerVector y_vec,
    int act_radius_signal, NumericMatrix signal_matrix, int grid_size)
{
  int n = x_vec.size();
  NumericVector results(n);

  for (int i = 0; i < n; i++)
  {
    results[i] = get_8n_avg_signal_cpp(
        x_vec[i], y_vec[i], act_radius_signal, signal_matrix, grid_size);
  }

  return results;
}

// ============================================================================
// MICROBE KILLING WITH ROS (VERY HOT LOOP)
// ============================================================================

// [[Rcpp::export]]
List kill_microbes_with_ros_cpp(
    NumericMatrix microbe_coords,
    NumericMatrix ROS,
    int act_radius_ROS,
    double th_ROS_microbe,
    int grid_size)
{
  int n_microbes = microbe_coords.nrow();

  if (n_microbes == 0)
  {
    NumericMatrix empty(0, 3);
    colnames(empty) = CharacterVector::create("x", "y", "id");
    return List::create(
        Named("surviving_microbes") = empty,
        Named("n_killed") = 0);
  }

  std::vector<int> survivors;
  survivors.reserve(n_microbes); // Pre-allocate for efficiency

  for (int i = 0; i < n_microbes; i++)
  {
    int x = static_cast<int>(microbe_coords(i, 0));
    int y = static_cast<int>(microbe_coords(i, 1));

    // Calculate average ROS using the fast function
    double avg_ros = get_8n_avg_signal_cpp(x, y, act_radius_ROS, ROS, grid_size);

    // Keep microbe if it survives
    if (avg_ros <= th_ROS_microbe)
    {
      survivors.push_back(i);
    }
  }

  // Build result matrix with survivors
  int n_survivors = survivors.size();
  NumericMatrix result(n_survivors, 3);

  for (int i = 0; i < n_survivors; i++)
  {
    int idx = survivors[i];
    result(i, 0) = microbe_coords(idx, 0);
    result(i, 1) = microbe_coords(idx, 1);
    result(i, 2) = microbe_coords(idx, 2);
  }

  colnames(result) = CharacterVector::create("x", "y", "id");

  return List::create(
      Named("surviving_microbes") = result,
      Named("n_killed") = n_microbes - n_survivors);
}

// ============================================================================
// PHAGOCYTE SIGNAL CALCULATIONS (HOT LOOP)
// ============================================================================

// Calculate DAMPs and SAMPs for all specified phagocytes
// [[Rcpp::export]]
List calculate_phagocyte_signals_cpp(
    IntegerVector phagocyte_indices,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    int act_radius_DAMPs,
    int act_radius_SAMPs,
    int act_radius_PAMPs,
    NumericMatrix DAMPs,
    NumericMatrix SAMPs,
    NumericMatrix PAMPs,
    int grid_size)
{
  int n = phagocyte_indices.size();
  NumericVector avg_DAMPs(n);
  NumericVector avg_SAMPs(n);
  NumericVector avg_PAMPs(n);
  // IntegerVector bacteria_counts(n);
  int bacteria_counts = 0;

  for (int i = 0; i < n; i++)
  {
    int idx = phagocyte_indices[i] - 1; // R to C++ indexing

    // Get position
    int x = phagocyte_x[idx];
    int y = phagocyte_y[idx];

    // Calculate average signals
    avg_DAMPs[i] = get_8n_avg_signal_cpp(x, y, act_radius_DAMPs, DAMPs, grid_size);
    avg_SAMPs[i] = get_8n_avg_signal_cpp(x, y, act_radius_SAMPs, SAMPs, grid_size);
    avg_PAMPs[i] = get_8n_avg_signal_cpp(x, y, act_radius_PAMPs, PAMPs, grid_size);

    // Count bacteria
    // int bacteria_count = 0;
    // for (int j = 0; j < phagocyte_bacteria_registry.ncol(); j++) {
    //   bacteria_count += phagocyte_bacteria_registry(idx, j);
    // }
    // bacteria_counts[i] = bacteria_count;
  }

  return List::create(
      Named("avg_DAMPs") = avg_DAMPs,
      Named("avg_SAMPs") = avg_SAMPs,
      Named("avg_PAMPs") = avg_PAMPs,
      Named("bacteria_counts") = bacteria_counts);
}

// ============================================================================
// ENGULFMENT PROCESSING (COMPLEX BOTTLENECK)
// ============================================================================

// Process engulfment for a single phagocyte
// Returns indices of microbes to remove
// [[Rcpp::export]]
IntegerVector find_overlapping_microbes_cpp(
    int px, int py,
    NumericMatrix microbe_coords)
{
  int n_microbes = microbe_coords.nrow();
  std::vector<int> overlapping;

  for (int i = 0; i < n_microbes; i++)
  {
    if (microbe_coords(i, 0) == px && microbe_coords(i, 1) == py)
    {
      overlapping.push_back(i + 1); // Convert to R indexing
    }
  }

  return wrap(overlapping);
}

// ============================================================================
// TREG VICINITY CALCULATIONS (HOT LOOP)
// ============================================================================

// Find tregs within vicinity of a phagocyte
// [[Rcpp::export]]
IntegerVector find_nearby_tregs_cpp(
    int px, int py,
    IntegerVector treg_x,
    IntegerVector treg_y,
    int treg_vicinity_effect)
{
  int n_tregs = treg_x.size();
  std::vector<int> nearby_indices;

  for (int i = 0; i < n_tregs; i++)
  {
    int dist_x = std::abs(treg_x[i] - px);
    int dist_y = std::abs(treg_y[i] - py);

    if (dist_x <= treg_vicinity_effect && dist_y <= treg_vicinity_effect)
    {
      nearby_indices.push_back(i + 1); // Convert to R indexing
    }
  }

  return wrap(nearby_indices);
}

// ============================================================================
// EPITHELIAL ROS CALCULATIONS (VECTORIZABLE)
// ============================================================================

// Calculate mean ROS for all epithelial cells at once
// [[Rcpp::export]]
NumericVector calculate_epithelial_ros_cpp(
    IntegerVector epithelium_x,
    int act_radius_ROS,
    NumericMatrix ROS,
    int grid_size)
{
  int n_cells = epithelium_x.size();
  NumericVector ros_means(n_cells);

  for (int i = 0; i < n_cells; i++)
  {
    int px = epithelium_x[i];
    int x0 = px - 1; // Convert to 0-based indexing

    // Calculate x coordinate range
    int x_start = std::max(0, x0 - act_radius_ROS);
    int x_end = std::min(grid_size - 1, x0 + act_radius_ROS);

    // ROS is at row 0 (epithelium is at y=1 in R, index 0 in C++)
    double sum = 0.0;
    int count = 0;

    for (int xi = x_start; xi <= x_end; xi++)
    {
      sum += ROS(0, xi); // Row 0 is the epithelium layer
      count++;
    }

    ros_means[i] = (count > 0) ? (sum / count) : 0.0;
  }

  return ros_means;
}

// ============================================================================
// MATRIX DIFFUSION (IF NEEDED)
// ============================================================================

// Optimized 8-neighbor diffusion
// [[Rcpp::export]]
NumericMatrix diffuse_matrix_cpp(
    NumericMatrix mat,
    double D,
    double max_cell_value,
    bool reflect_top = false)
{
  int nr = mat.nrow();
  int nc = mat.ncol();

  // Create padded matrix
  NumericMatrix padded(nr + 2, nc + 2);
  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      padded(i + 1, j + 1) = mat(i, j);
    }
  }

  // REFLECTIVE boundary at TOP (y=0, epithelium boundary)
  if (reflect_top)
  {
    for (int j = 0; j < nc; j++)
    {
      padded(0, j + 1) = mat(0, j); // Reflect epithelium
    }
  }

  // Calculate 8-neighbor Laplacian
  NumericMatrix result(nr, nc);

  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      double laplacian =
          padded(i, j) +         // top-left
          padded(i, j + 1) +     // top
          padded(i, j + 2) +     // top-right
          padded(i + 1, j) +     // left
          padded(i + 1, j + 2) + // right
          padded(i + 2, j) +     // bottom-left
          padded(i + 2, j + 1) + // bottom
          padded(i + 2, j + 2) - // bottom-right
          8.0 * mat(i, j);       // center

      double new_val = mat(i, j) + D * laplacian;
      result(i, j) = std::min(max_cell_value, new_val);
    }
  }

  return result;
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

// Fast shift and insert for bacteria registry
// [[Rcpp::export]]
NumericVector shift_insert_fast_cpp(
    NumericVector vec,
    NumericVector insert_vals)
{
  int n_insert = insert_vals.size();
  int n_vec = vec.size();

  NumericVector result(n_vec);

  if (n_insert >= n_vec)
  {
    // If inserting more than vec size, just take first n_vec elements
    for (int i = 0; i < n_vec; i++)
    {
      result[i] = insert_vals[i];
    }
  }
  else
  {
    // Insert at beginning, shift rest
    for (int i = 0; i < n_insert; i++)
    {
      result[i] = insert_vals[i];
    }
    for (int i = n_insert; i < n_vec; i++)
    {
      result[i] = vec[i - n_insert];
    }
  }

  return result;
}

// Optimized coordinate movement helper
// [[Rcpp::export]]
IntegerVector iszero_coordinates_cpp(IntegerVector x)
{
  int n = x.size();
  IntegerVector y(n);

  // First pass: sample for all positions
  NumericVector rand_vals = runif(n, -1, 1);

  for (int i = 0; i < n; i++)
  {
    if (x[i] == 0)
    {
      // For zero values, choose -1 or 1
      y[i] = (rand_vals[i] < 0) ? -1 : 1;
    }
    else
    {
      // For non-zero, sample from {-1, 0, 1}
      double r = runif(1)[0];
      if (r < 0.333)
      {
        y[i] = -1;
      }
      else if (r < 0.667)
      {
        y[i] = 0;
      }
      else
      {
        y[i] = 1;
      }
    }
  }

  return y;
}

// ============================================================================
// CHEMOTAXIS MOVEMENT (MAJOR OPTIMIZATION)
// ============================================================================

// Move all cells at once using chemotaxis gradient
// [[Rcpp::export]]
List move_cells_chemotaxis_cpp(
    IntegerVector cell_x,
    IntegerVector cell_y,
    NumericMatrix density_matrix,
    int grid_size)
{
  int n_cells = cell_x.size();
  IntegerVector new_x(n_cells);
  IntegerVector new_y(n_cells);

  for (int i = 0; i < n_cells; i++)
  {
    int x = cell_x[i];
    int y = cell_y[i];

    // Calculate bounds (R indexing)
    int x_min = std::max(1, x - 1);
    int x_max = std::min(grid_size, x + 1);
    int y_min = std::max(1, y - 1);
    int y_max = std::min(grid_size, y + 1);

    // Collect neighbors
    std::vector<int> neighbors_x;
    std::vector<int> neighbors_y;
    std::vector<double> neighbor_densities;

    for (int yi = y_min; yi <= y_max; yi++)
    {
      for (int xi = x_min; xi <= x_max; xi++)
      {
        neighbors_x.push_back(xi);
        neighbors_y.push_back(yi);
        // Convert to 0-based for matrix access
        neighbor_densities.push_back(density_matrix(yi - 1, xi - 1));
      }
    }

    // Calculate probabilities
    double total = 0.0;
    for (double d : neighbor_densities)
    {
      total += d;
    }

    // Sample destination
    int chosen_idx = 0;
    if (total > 0)
    {
      // Weighted sampling
      NumericVector probs(neighbor_densities.size());
      for (size_t j = 0; j < neighbor_densities.size(); j++)
      {
        probs[j] = neighbor_densities[j] / total;
      }

      double u = runif(1)[0];
      double cumsum = 0.0;
      for (size_t j = 0; j < probs.size(); j++)
      {
        cumsum += probs[j];
        if (u <= cumsum)
        {
          chosen_idx = j;
          break;
        }
      }
    }
    else
    {
      // Uniform random if all densities are 0
      chosen_idx = floor(runif(1)[0] * neighbors_x.size());
    }

    new_x[i] = neighbors_x[chosen_idx];
    new_y[i] = neighbors_y[chosen_idx];
  }

  return List::create(
      Named("x") = new_x,
      Named("y") = new_y);
}

// ============================================================================
// ENGULFMENT PROCESS (MAJOR OPTIMIZATION)
// ============================================================================

// Process all engulfments at once for all phagocytes
// [[Rcpp::export]]
List process_engulfment_batch_cpp(
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    NumericVector phagocyte_activity_engulf,
    IntegerVector phagocyte_phenotype,
    NumericMatrix pathogen_coords,
    NumericMatrix commensal_coords,
    NumericMatrix phagocyte_bacteria_registry,
    int cc_phagocyte)
{
  int n_phagocytes = phagocyte_x.size();
  int n_pathogens = pathogen_coords.nrow();
  int n_commensals = commensal_coords.nrow();

  // Track which microbes were engulfed
  std::vector<bool> pathogen_alive(n_pathogens, true);
  std::vector<bool> commensal_alive(n_commensals, true);

  // Track phagocyte engulfment success
  LogicalVector phagocytes_that_engulfed(n_phagocytes, false);

  // Updated registry
  NumericMatrix updated_registry = clone(phagocyte_bacteria_registry);

  // Kill counters [M0, M1, M2]
  IntegerVector pathogens_killed(3, 0);
  IntegerVector commensals_killed(3, 0);

  // Process each phagocyte
  for (int i = 0; i < n_phagocytes; i++)
  {
    int px = phagocyte_x[i];
    int py = phagocyte_y[i];
    double engulf_prob = phagocyte_activity_engulf[i];
    int phenotype = phagocyte_phenotype[i];

    std::vector<int> engulfed_values;

    // Check pathogen overlaps
    for (int p = 0; p < n_pathogens; p++)
    {
      if (pathogen_alive[p] &&
          pathogen_coords(p, 0) == px &&
          pathogen_coords(p, 1) == py)
      {
        // Try to engulf
        if (runif(1)[0] < engulf_prob)
        {
          pathogen_alive[p] = false;
          engulfed_values.push_back(-1); // -1 for pathogen
          pathogens_killed[phenotype]++;
        }
      }
    }

    // Check commensal overlaps
    for (int c = 0; c < n_commensals; c++)
    {
      if (commensal_alive[c] &&
          commensal_coords(c, 0) == px &&
          commensal_coords(c, 1) == py)
      {
        // Try to engulf
        if (runif(1)[0] < engulf_prob)
        {
          commensal_alive[c] = false;
          engulfed_values.push_back(1); // +1 for commensal
          commensals_killed[phenotype]++;
        }
      }
    }

    // Update registry if engulfment occurred
    if (engulfed_values.size() > 0)
    {
      phagocytes_that_engulfed[i] = true;

      // Shift and insert
      int n_insert = engulfed_values.size();
      if (n_insert >= cc_phagocyte)
      {
        for (int j = 0; j < cc_phagocyte; j++)
        {
          updated_registry(i, j) = engulfed_values[j];
        }
      }
      else
      {
        // Insert at beginning, shift rest
        for (int j = 0; j < n_insert; j++)
        {
          updated_registry(i, j) = engulfed_values[j];
        }
        for (int j = n_insert; j < cc_phagocyte; j++)
        {
          updated_registry(i, j) = phagocyte_bacteria_registry(i, j - n_insert);
        }
      }
    }
  }

  // Build surviving microbe matrices
  std::vector<int> surviving_pathogen_idx;
  std::vector<int> surviving_commensal_idx;

  for (int p = 0; p < n_pathogens; p++)
  {
    if (pathogen_alive[p])
      surviving_pathogen_idx.push_back(p);
  }
  for (int c = 0; c < n_commensals; c++)
  {
    if (commensal_alive[c])
      surviving_commensal_idx.push_back(c);
  }

  NumericMatrix surviving_pathogens(surviving_pathogen_idx.size(), 3);
  NumericMatrix surviving_commensals(surviving_commensal_idx.size(), 3);

  for (size_t i = 0; i < surviving_pathogen_idx.size(); i++)
  {
    int idx = surviving_pathogen_idx[i];
    for (int j = 0; j < 3; j++)
    {
      surviving_pathogens(i, j) = pathogen_coords(idx, j);
    }
  }
  colnames(surviving_pathogens) = CharacterVector::create("x", "y", "id");

  for (size_t i = 0; i < surviving_commensal_idx.size(); i++)
  {
    int idx = surviving_commensal_idx[i];
    for (int j = 0; j < 3; j++)
    {
      surviving_commensals(i, j) = commensal_coords(idx, j);
    }
  }
  colnames(surviving_commensals) = CharacterVector::create("x", "y", "id");

  return List::create(
      Named("pathogen_coords") = surviving_pathogens,
      Named("commensal_coords") = surviving_commensals,
      Named("phagocyte_bacteria_registry") = updated_registry,
      Named("phagocytes_that_engulfed") = phagocytes_that_engulfed,
      Named("pathogens_killed") = pathogens_killed,
      Named("commensals_killed") = commensals_killed);
}

// ============================================================================
// PAMPs UPDATE (HOT LOOP OPTIMIZATION)
// ============================================================================

// Update PAMPs matrix from pathogen coordinates
// [[Rcpp::export]]
NumericMatrix update_PAMPs_cpp(
    NumericMatrix PAMPs,
    NumericMatrix pathogen_coords,
    double add_PAMPs,
    double k_in_log,
    double x0_in_log,
    double c_in_log)
{
  NumericMatrix result = clone(PAMPs);

  if (pathogen_coords.nrow() == 0)
  {
    return result;
  }

  // Count pathogens at each location
  std::map<std::pair<int, int>, int> pathogen_counts;

  for (int i = 0; i < pathogen_coords.nrow(); i++)
  {
    int x = static_cast<int>(pathogen_coords(i, 0));
    int y = static_cast<int>(pathogen_coords(i, 1));
    pathogen_counts[std::make_pair(x, y)]++;
  }

  // Add PAMPs for each location
  for (auto &kv : pathogen_counts)
  {
    int x = kv.first.first;
    int y = kv.first.second;
    int count = kv.second;

    // logistic_scaled_0_to_5_quantized
    double logistic_val = c_in_log / (1.0 + exp(-k_in_log * (count - x0_in_log)));
    double pamps_to_add = add_PAMPs * logistic_val;

    // Convert to 0-based indexing
    result(y - 1, x - 1) += pamps_to_add;
  }

  return result;
}

// ============================================================================
// TREG ACTIVATION BATCH (OPTIMIZATION)
// ============================================================================

// Process treg activation for all M1/M2 phagocytes at once
// [[Rcpp::export]]
List activate_tregs_batch_cpp(
    IntegerVector M_phagocyte_indices,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    IntegerVector phagocyte_pathogens_engulfed,
    IntegerVector phagocyte_commensals_engulfed,
    IntegerVector treg_x,
    IntegerVector treg_y,
    IntegerVector treg_phenotype,
    NumericVector treg_activity_SAMPs_binary,
    IntegerVector treg_active_age,
    int treg_vicinity_effect,
    double treg_discrimination_efficiency)
{
  int n_tregs = treg_x.size();

  // Clone outputs
  IntegerVector new_treg_phenotype = clone(treg_phenotype);
  NumericVector new_treg_activity = clone(treg_activity_SAMPs_binary);
  IntegerVector new_treg_age = clone(treg_active_age);

  // Process each activated phagocyte
  for (int idx = 0; idx < M_phagocyte_indices.size(); idx++)
  {
    int i = M_phagocyte_indices[idx] - 1; // R to C++ indexing
    int px = phagocyte_x[i];
    int py = phagocyte_y[i];

    int num_pat = phagocyte_pathogens_engulfed[i];
    int num_com = phagocyte_commensals_engulfed[i];

    if ((num_pat + num_com) == 0)
      continue;

    // Find nearby tregs
    std::vector<int> nearby_tregs;
    for (int t = 0; t < n_tregs; t++)
    {
      int dist_x = std::abs(treg_x[t] - px);
      int dist_y = std::abs(treg_y[t] - py);

      if (dist_x <= treg_vicinity_effect && dist_y <= treg_vicinity_effect)
      {
        nearby_tregs.push_back(t);
      }
    }

    if (nearby_tregs.size() == 0)
      continue;

    // Calculate probability of presenting commensal
    double rat_com_pat_real = (double)num_com / (num_com + num_pat);

    // Step 1: Which antigen is presented?
    bool commensal_presented = runif(1)[0] < rat_com_pat_real;

    // Step 2: Does Treg identify it correctly?
    bool treg_identifies_as_commensal;
    if (commensal_presented)
    {
      treg_identifies_as_commensal = runif(1)[0] < treg_discrimination_efficiency;
    }
    else
    {
      treg_identifies_as_commensal = runif(1)[0] < (1.0 - treg_discrimination_efficiency);
    }

    // Step 3: Activate if identified as commensal
    if (treg_identifies_as_commensal)
    {
      for (int t_idx : nearby_tregs)
      {
        new_treg_phenotype[t_idx] = 1;
        new_treg_activity[t_idx] = 1;
        new_treg_age[t_idx] = 1;
      }
    }
  }

  return List::create(
      Named("treg_phenotype") = new_treg_phenotype,
      Named("treg_activity_SAMPs_binary") = new_treg_activity,
      Named("treg_active_age") = new_treg_age);
}

// ============================================================================
// PHAGOCYTE PHENOTYPE TRANSITIONS (BATCH OPTIMIZATION)
// ============================================================================

// Process M0 -> M1/M2 transitions
// [[Rcpp::export]]
List update_M0_phenotypes_cpp(
    IntegerVector M0_indices,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    NumericMatrix DAMPs,
    NumericMatrix SAMPs,
    NumericMatrix PAMPs,
    int act_radius_DAMPs,
    int act_radius_SAMPs,
    int act_radius_PAMPs,
    int grid_size,
    double activation_threshold_SAMPs,
    double activation_threshold_danger,
    double activity_ROS_M1_baseline,
    double activity_engulf_M1_baseline)
{
  int n = M0_indices.size();

  IntegerVector new_phenotype(n, 0);
  IntegerVector new_active_age(n, 0);
  NumericVector new_activity_ROS(n);
  NumericVector new_activity_engulf(n);
  LogicalVector activated(n, false);

  for (int idx = 0; idx < n; idx++)
  {
    int i = M0_indices[idx] - 1; // R to C++ indexing
    int x = phagocyte_x[i];
    int y = phagocyte_y[i];

    // Calculate signals
    double avg_DAMPs = get_8n_avg_signal_cpp(x, y, act_radius_DAMPs, DAMPs, grid_size);
    double avg_SAMPs = get_8n_avg_signal_cpp(x, y, act_radius_SAMPs, SAMPs, grid_size);
    double avg_PAMPs = get_8n_avg_signal_cpp(x, y, act_radius_PAMPs, PAMPs, grid_size);

    double danger_signal = avg_DAMPs + avg_PAMPs;

    double SAMPS_diff = std::max(0.0, avg_SAMPs - activation_threshold_SAMPs) / activation_threshold_SAMPs;
    double DANGER_diff = std::max(0.0, danger_signal - activation_threshold_danger) / activation_threshold_danger;

    if (SAMPS_diff > 0 || DANGER_diff > 0)
    {
      if (DANGER_diff > SAMPS_diff)
      {
        // M0 -> M1
        new_phenotype[idx] = 1;
        new_active_age[idx] = 1;
        new_activity_ROS[idx] = activity_ROS_M1_baseline;
        new_activity_engulf[idx] = activity_engulf_M1_baseline;
        activated[idx] = true;
      }
    }
  }

  return List::create(
      Named("phenotype") = new_phenotype,
      Named("active_age") = new_active_age,
      Named("activity_ROS") = new_activity_ROS,
      Named("activity_engulf") = new_activity_engulf,
      Named("activated") = activated);
}

// Process M1/M2 phenotype transitions and deactivation
// [[Rcpp::export]]
List update_active_phenotypes_cpp(
    IntegerVector active_indices,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    IntegerVector phagocyte_phenotype,
    IntegerVector phagocyte_active_age,
    NumericMatrix DAMPs,
    NumericMatrix SAMPs,
    NumericMatrix PAMPs,
    int act_radius_DAMPs,
    int act_radius_SAMPs,
    int act_radius_PAMPs,
    int grid_size,
    int active_age_limit,
    double activation_threshold_SAMPs,
    double activation_threshold_danger,
    double activity_ROS_M0_baseline,
    double activity_engulf_M0_baseline,
    double activity_ROS_M1_baseline,
    double activity_engulf_M1_baseline,
    double activity_ROS_M2_baseline,
    double activity_engulf_M2_baseline,
    bool only_suppress)
{
  int n = active_indices.size();

  IntegerVector new_phenotype(n);
  IntegerVector new_active_age(n);
  NumericVector new_activity_ROS(n);
  NumericVector new_activity_engulf(n);
  LogicalVector changed(n, false);

  for (int idx = 0; idx < n; idx++)
  {
    int i = active_indices[idx] - 1; // R to C++ indexing
    int x = phagocyte_x[i];
    int y = phagocyte_y[i];
    int current_phenotype = phagocyte_phenotype[i];
    int current_age = phagocyte_active_age[i];

    // Default: no change
    new_phenotype[idx] = current_phenotype;
    new_active_age[idx] = current_age;
    new_activity_ROS[idx] = (current_phenotype == 1) ? activity_ROS_M1_baseline : activity_ROS_M2_baseline;
    new_activity_engulf[idx] = (current_phenotype == 1) ? activity_engulf_M1_baseline : activity_engulf_M2_baseline;

    // Check if old enough for transition
    bool old_enough = current_age >= active_age_limit;

    if (!old_enough && only_suppress)
    {
      // Young cells can only suppress (M1->M2 or M2->M1)
      double avg_DAMPs = get_8n_avg_signal_cpp(x, y, act_radius_DAMPs, DAMPs, grid_size);
      double avg_SAMPs = get_8n_avg_signal_cpp(x, y, act_radius_SAMPs, SAMPs, grid_size);
      double avg_PAMPs = get_8n_avg_signal_cpp(x, y, act_radius_PAMPs, PAMPs, grid_size);

      double danger_signal = avg_DAMPs + avg_PAMPs;

      double SAMPS_diff = std::max(0.0, avg_SAMPs - activation_threshold_SAMPs) / activation_threshold_SAMPs;
      double DANGER_diff = std::max(0.0, danger_signal - activation_threshold_danger) / activation_threshold_danger;

      if (current_phenotype == 1 && SAMPS_diff > DANGER_diff)
      {
        // M1 -> M2
        new_phenotype[idx] = 2;
        new_active_age[idx] = 1;
        new_activity_ROS[idx] = activity_ROS_M2_baseline;
        new_activity_engulf[idx] = activity_engulf_M2_baseline;
        changed[idx] = true;
      }
      else if (current_phenotype == 2 && SAMPS_diff < DANGER_diff)
      {
        // M2 -> M1
        new_phenotype[idx] = 1;
        new_active_age[idx] = 1;
        new_activity_ROS[idx] = activity_ROS_M1_baseline;
        new_activity_engulf[idx] = activity_engulf_M1_baseline;
        changed[idx] = true;
      }
    }
    else if (old_enough && !only_suppress)
    {
      // Old cells can deactivate or transition
      double avg_DAMPs = get_8n_avg_signal_cpp(x, y, act_radius_DAMPs, DAMPs, grid_size);
      double avg_SAMPs = get_8n_avg_signal_cpp(x, y, act_radius_SAMPs, SAMPs, grid_size);
      double avg_PAMPs = get_8n_avg_signal_cpp(x, y, act_radius_PAMPs, PAMPs, grid_size);

      double danger_signal = avg_DAMPs + avg_PAMPs;

      double SAMPS_diff = std::max(0.0, avg_SAMPs - activation_threshold_SAMPs) / activation_threshold_SAMPs;
      double DANGER_diff = std::max(0.0, danger_signal - activation_threshold_danger) / activation_threshold_danger;

      if (current_phenotype == 1)
      {
        if (SAMPS_diff == 0 && DANGER_diff == 0)
        {
          // M1 -> M0
          new_phenotype[idx] = 0;
          new_active_age[idx] = 0;
          new_activity_ROS[idx] = activity_ROS_M0_baseline;
          new_activity_engulf[idx] = activity_engulf_M0_baseline;
          changed[idx] = true;
        }
        else if (SAMPS_diff > DANGER_diff)
        {
          // M1 -> M2
          new_phenotype[idx] = 2;
          new_active_age[idx] = 1;
          new_activity_ROS[idx] = activity_ROS_M2_baseline;
          new_activity_engulf[idx] = activity_engulf_M2_baseline;
          changed[idx] = true;
        }
      }
      else if (current_phenotype == 2)
      {
        if (SAMPS_diff == 0 && DANGER_diff == 0)
        {
          // M2 -> M0
          new_phenotype[idx] = 0;
          new_active_age[idx] = 0;
          new_activity_ROS[idx] = activity_ROS_M0_baseline;
          new_activity_engulf[idx] = activity_engulf_M0_baseline;
          changed[idx] = true;
        }
        else if (SAMPS_diff < DANGER_diff)
        {
          // M2 -> M1
          new_phenotype[idx] = 1;
          new_active_age[idx] = 1;
          new_activity_ROS[idx] = activity_ROS_M1_baseline;
          new_activity_engulf[idx] = activity_engulf_M1_baseline;
          changed[idx] = true;
        }
      }
    }
  }

  return List::create(
      Named("phenotype") = new_phenotype,
      Named("active_age") = new_active_age,
      Named("activity_ROS") = new_activity_ROS,
      Named("activity_engulf") = new_activity_engulf,
      Named("changed") = changed);
}

// ============================================================================
// REGISTRY SHIFT VECTORIZED
// ============================================================================

// Shift all registries at once for phagocytes that didn't engulf
// [[Rcpp::export]]
NumericMatrix shift_registry_batch_cpp(
    NumericMatrix registry,
    LogicalVector should_shift)
{
  int n_rows = registry.nrow();
  int n_cols = registry.ncol();

  NumericMatrix result = clone(registry);

  for (int i = 0; i < n_rows; i++)
  {
    if (should_shift[i])
    {
      // Shift left (older bacteria move to the right)
      for (int j = n_cols - 1; j > 0; j--)
      {
        result(i, j) = registry(i, j - 1);
      }
      result(i, 0) = 0; // Empty slot at the front
    }
  }

  return result;
}

// ============================================================================
// BATCH OPERATIONS FOR MAXIMUM EFFICIENCY
// ============================================================================

// Update SAMPs matrix with multiple active tregs at once
// [[Rcpp::export]]
NumericMatrix update_SAMPs_batch_cpp(
    NumericMatrix SAMPs,
    IntegerVector active_tregs,
    IntegerVector treg_x,
    IntegerVector treg_y,
    NumericVector treg_activity_SAMPs_binary,
    double add_SAMPs,
    int allow_tregs)
{
  NumericMatrix result = clone(SAMPs);
  int n_active = active_tregs.size();

  for (int i = 0; i < n_active; i++)
  {
    int idx = active_tregs[i] - 1; // R to C++ indexing
    int x = treg_x[idx] - 1;
    int y = treg_y[idx] - 1;

    result(y, x) += treg_activity_SAMPs_binary[idx] * add_SAMPs * allow_tregs;
  }

  return result;
}
