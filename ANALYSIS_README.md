# Generic ABM Region Analysis Script

## Overview

The `DLL_datazanalyse_abm_regions_generic.R` script consolidates three separate analysis scripts into one flexible, parameterized script that can handle any pairwise comparison between treatment conditions.

## Original Scripts Consolidated

1. `DLL_datazanalyse_abm_regions_macspec_vs_treg_off.R`
2. `DLL_datazanalyse_abm_regions_tregs_on_va_macspec.R`
3. `DLL_datazanalyse_abm_regions_tregs_on_vs_tregs_off.R`

## Quick Start

### Option 1: Run All Original Analyses

```r
source("./run_all_comparisons.R")
```

This will reproduce all three original analyses with a single command.

### Option 2: Run Individual Comparisons

Set configuration parameters at the top of the generic script or in your R session:

```r
# Set parameters
group1 <- "tregs_on"
group2 <- "tregs_off"
difference_direction <- "group1_minus_group2"
filter_control <- 1
jsd_th <- 0.3
tol_in <- 18.75
output_prefix <- "tregs_on_vs_tregs_off"
check_randomization_logic <- TRUE
show_intermediate_picks <- FALSE

# Run analysis
source("./DLL_datazanalyse_abm_regions_generic.R")
```

## Configuration Parameters

### Required Parameters

| Parameter | Description | Example Values |
|-----------|-------------|----------------|
| `group1` | First treatment group | `"ctrl"`, `"tregs_off"`, `"tregs_on"`, `"tregs_rnd"`, `"macspec"` |
| `group2` | Second treatment group | `"ctrl"`, `"tregs_off"`, `"tregs_on"`, `"tregs_rnd"`, `"macspec"` |
| `difference_direction` | How to calculate difference | `"group1_minus_group2"` or `"group2_minus_group1"` |
| `output_prefix` | Prefix for output files | `"tregs_on_vs_tregs_off"` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filter_control` | `0` | Apply control filtering: `1` = yes, `0` = no |
| `jsd_th` | `0.3` | Cohen's d threshold for significance |
| `tol_in` | `18.75` | Tolerance threshold for differences |
| `check_randomization_logic` | `FALSE` | Check randomization effects (for tregs_on vs tregs_off) |
| `show_intermediate_picks` | `FALSE` | Show intermediate filtering statistics |

## Available Treatment Groups

The script supports comparisons between:

- `ctrl` - Control condition
- `tregs_off` - Tregs turned off
- `tregs_on` - Tregs turned on
- `tregs_rnd` - Tregs randomized
- `macspec` - Macrophage specific condition

## Output Files

For each analysis, the script generates:

1. **RDS file**: `df_comparisons_plot_<output_prefix>.rds`
   - Contains the full comparison data frame
   - Can be loaded for further analysis

2. **PNG file**: `ABM_JSD_<OUTPUT_PREFIX>.png`
   - Visualization of the comparison
   - Shows Cohen's d vs difference plot with regions

## Examples

### Example 1: Compare tregs_on vs ctrl

```r
group1 <- "tregs_on"
group2 <- "ctrl"
difference_direction <- "group1_minus_group2"
filter_control <- 0
output_prefix <- "tregs_on_vs_ctrl"
check_randomization_logic <- FALSE
show_intermediate_picks <- FALSE

source("./DLL_datazanalyse_abm_regions_generic.R")
```

### Example 2: Compare macspec vs tregs_on with control filter

```r
group1 <- "macspec"
group2 <- "tregs_on"
difference_direction <- "group1_minus_group2"
filter_control <- 1
output_prefix <- "macspec_vs_tregs_on"
check_randomization_logic <- FALSE
show_intermediate_picks <- TRUE

source("./DLL_datazanalyse_abm_regions_generic.R")
```

### Example 3: Reverse comparison (group2 - group1)

```r
group1 <- "tregs_off"
group2 <- "macspec"
difference_direction <- "group2_minus_group1"
filter_control <- 0
output_prefix <- "macspec_vs_tregs_off_reversed"
check_randomization_logic <- FALSE
show_intermediate_picks <- FALSE

source("./DLL_datazanalyse_abm_regions_generic.R")
```

## How It Works

### 1. Data Loading and Filtering
- Loads parameters and simulation results
- Filters for complete replication sets
- Filters based on `ss_start` threshold

### 2. Control Filtering (Optional)
- If `filter_control = 1`, filters based on control pathogenic comparison
- Keeps only parameter sets where control differs significantly from test

### 3. Difference Calculation
- Dynamically constructs variable names based on `group1` and `group2`
- Calculates differences for both pathogenic and sterile injuries
- Direction determined by `difference_direction` parameter

### 4. Cohen's d Analysis
- Automatically finds correct Cohen's d columns
- Handles different naming conventions
- Separates analysis for pathogenic and sterile conditions

### 5. Conflict Detection
- Identifies parameter sets where effect differs between injury types
- Reports cases where treatment is better/worse in different conditions

### 6. Output Generation
- Saves data frame and generates visualization
- Provides summary statistics

## Understanding the Output

### Summary Statistics

The script reports:

- **Pathogenic injury**: % better, % worse
- **Sterile injury**: % better, % worse
- **Both injuries**: % better in both conditions
- **Conflicts**: Cases where direction differs between injury types

### Plot Interpretation

The generated plot shows:

- **X-axis**: Cohen's d (effect size)
- **Y-axis**: Mean difference in score
- **Regions**: Different colors indicate beneficial/detrimental/neutral effects
- **Points**: Each parameter set for pathogenic and sterile conditions

## Troubleshooting

### Column Not Found Error

If you see an error like "Could not find Cohen's d columns", check:

1. Your group names match exactly what's in the data
2. The comparison exists in the original simulation data
3. Try reversing `group1` and `group2`

### No Data After Filtering

If no data remains after filtering:

- Reduce `jsd_th` threshold
- Set `filter_control = 0`
- Check `tol_in` value isn't too restrictive

## Extending the Script

To add new comparisons:

1. Ensure mean score columns exist in data: `mean_<group>_<injury>_<metric>`
2. Ensure Cohen's d columns exist: `d_<group1>_vs_<group2>_<injury>_<metric>`
3. Set appropriate parameters and source the generic script

## Dependencies

Required R packages:
- dplyr
- tidyr
- ggplot2
- purrr
- readr
- stringr
- zoo
- scales
- ggrepel

Required external scripts:
- `./MISC/PLOT_FUNCTIONS_ABM.R`
- `./MISC/DATA_READ_FUNCTIONS.R`
- `./MISC/REGIONS.R`

Required data files:
- `./lhs_parameters_della.csv`
- `./data_cpp_read_abm.rds`
