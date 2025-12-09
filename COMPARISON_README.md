# Generic Comparison Script for ABM Regions Analysis

## Overview

This system provides a generic, reusable framework for comparing different experimental groups in the ABM regions analysis. Instead of maintaining multiple similar scripts, you can now use a single generic script (`DLL_datazanalyse_abm_regions_generic.R`) configured through parameters.

## Files

- **`DLL_datazanalyse_abm_regions_generic.R`**: The generic comparison script that performs all analyses
- **`run_comparisons.R`**: Configuration file with predefined comparison functions
- **Original scripts** (now superseded by the generic approach):
  - `DLL_datazanalyse_abm_regions_tregs_on_vs_tregs_off.R`
  - `DLL_datazanalyse_abm_regions_tregs_on_vs_macspec.R`
  - `DLL_datazanalyse_abm_regions_macspec_vs_treg_off.R`

## Quick Start

### Running Predefined Comparisons

```r
# Load the configuration
source('./run_comparisons.R')

# Run individual comparisons
run_comparison_tregs_on_vs_off()
run_comparison_tregs_on_vs_macspec()
run_comparison_macspec_vs_tregs_off()

# Or run all comparisons at once
run_all_comparisons()
```

### Creating Custom Comparisons

```r
source('./run_comparisons.R')

# Custom comparison between any two groups
run_comparison(
  comparison_name = "my_custom_comparison",
  group1 = "tregs_on",
  group2 = "ctrl",
  jsd_th = 0.3,
  tol_in = 125 * 0.2,
  filter_control = 0,
  perform_randomness_check = FALSE,
  perform_macspec_analysis = FALSE
)
```

## Parameters

### Required Parameters

- **`comparison_name`**: Identifier for the comparison (used in output filenames)
- **`group1`**: Name of the first group (e.g., "tregs_on", "macspec", "ctrl")
- **`group2`**: Name of the second group (e.g., "tregs_off", "macspec")
- **`jsd_th`**: Jensen-Shannon distance threshold for significance
- **`filter_control`**: Whether to filter based on control comparison (0 or 1)

### Optional Parameters

- **`tol_in`**: Tolerance threshold for difference (set to `NA` for auto-calculation)
- **`perform_randomness_check`**: Include analysis of randomized tregs (default: `FALSE`)
- **`perform_macspec_analysis`**: Include macspec-specific analysis (default: `FALSE`)

## Available Groups

The following groups are available for comparison:
- `ctrl`: Control group
- `tregs_off`: Tregs off condition
- `tregs_on`: Tregs on condition
- `tregs_rnd`: Randomized tregs
- `macspec`: Macrophage-specific condition

## Outputs

For each comparison, the script generates:

1. **RDS file**: `df_comparisons_plot_<comparison_name>.rds`
   - Contains processed comparison data
   - Includes Cohen's d values, significance flags, and region classifications

2. **PNG plot**: `ABM_JSD_<COMPARISON_NAME>.png`
   - Visualization of the comparison across parameter space
   - Shows pathogenic vs sterile injury types
   - Color-coded by effect magnitude and direction

3. **Console output**: Summary statistics including:
   - Number of parameter sets analyzed
   - Percentage better/worse for pathogenic injury
   - Percentage better/worse for sterile injury
   - Percentage better for both injury types
   - Conflicting cases (if any)

## How It Works

### Data Processing Pipeline

1. **Data Loading**: Loads parameter sets and simulation results
2. **Filtering**: Filters for complete replicates and appropriate starting states
3. **Comparison Setup**: Selects relevant mean scores for both groups
4. **Control Filtering** (optional): Filters based on control group comparison
5. **Difference Calculation**: Computes differences between groups
6. **Statistical Testing**: Merges with Cohen's d effect sizes
7. **Significance Classification**: Determines which cases show significant effects
8. **Conflict Analysis**: Identifies cases with opposite effects in sterile vs pathogenic
9. **Visualization**: Generates publication-ready plots
10. **Reporting**: Outputs summary statistics

### Cohen's D Column Naming

The script automatically handles different naming conventions for Cohen's d columns:
- Tries `d_<group1>_vs_<group2>_<injury>_<metric>`
- Falls back to `d_<group2>_vs_<group1>_<injury>_<metric>` if needed

## Examples

### Example 1: Replicate Original Analysis

```r
source('./run_comparisons.R')

# This replicates DLL_datazanalyse_abm_regions_tregs_on_vs_tregs_off.R
run_comparison_tregs_on_vs_off()
```

### Example 2: Compare Control vs Macspec

```r
source('./run_comparisons.R')

run_comparison(
  comparison_name = "ctrl_vs_macspec",
  group1 = "ctrl",
  group2 = "macspec",
  jsd_th = 0.3,
  tol_in = 125 * 0.15,
  filter_control = 0
)
```

### Example 3: Batch Processing Multiple Comparisons

```r
source('./run_comparisons.R')

# Define comparisons of interest
comparisons = list(
  list(name = "tregs_on_vs_off", g1 = "tregs_on", g2 = "tregs_off"),
  list(name = "tregs_on_vs_rnd", g1 = "tregs_on", g2 = "tregs_rnd"),
  list(name = "ctrl_vs_tregs_on", g1 = "ctrl", g2 = "tregs_on")
)

# Run all comparisons
for (comp in comparisons) {
  run_comparison(
    comparison_name = comp$name,
    group1 = comp$g1,
    group2 = comp$g2,
    jsd_th = 0.3,
    tol_in = 125 * 0.2,
    filter_control = 0
  )
}
```

## Benefits Over Original Scripts

1. **Single Source of Truth**: All comparison logic in one place
2. **Easier Maintenance**: Bug fixes and improvements apply to all comparisons
3. **Consistent Results**: Same methodology across all comparisons
4. **Flexible**: Easy to add new comparisons without duplicating code
5. **Well-Documented**: Clear parameter definitions and expected outputs
6. **Reproducible**: Explicit parameter settings for each comparison

## Troubleshooting

### Missing Cohen's D Column

If you get an error about missing Cohen's d columns, check that:
1. The group names match exactly those in the data
2. The Cohen's d columns exist in `df_results`
3. The naming convention matches `d_<group1>_vs_<group2>_<injury>_e`

### Auto-calculated Tolerance Too Large/Small

If the auto-calculated `tol_in` seems inappropriate:
1. Set it manually based on domain knowledge
2. Check the distribution of differences in non-significant cases
3. Consider using a fraction of the maximum possible score (e.g., `125 * 0.2`)

### No Significant Results

If no parameter sets show significant effects:
1. Try lowering `jsd_th` (Cohen's d threshold)
2. Try lowering `tol_in` (difference threshold)
3. Check if `filter_control` is too restrictive

## Migration from Original Scripts

If you were using the original scripts, simply replace:

```r
# Old approach
source('./DLL_datazanalyse_abm_regions_tregs_on_vs_tregs_off.R')
```

With:

```r
# New approach
source('./run_comparisons.R')
run_comparison_tregs_on_vs_off()
```

The results should be identical (except for potential minor differences due to the fixed column naming logic).
