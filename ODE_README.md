# ODE Version of Immune Response Model

This directory contains an Ordinary Differential Equation (ODE) version of the agent-based immune response model, designed to capture population-level dynamics without spatial structure or individual-level stochasticity.

## Overview

The ODE system translates the spatial agent-based model (ABM) into a well-mixed, deterministic system with 10 state variables:

### State Variables (Compartments)
- **M0, M1, M2**: Macrophage populations (M0 = naive, M1 = pro-inflammatory, M2 = anti-inflammatory)
- **P**: Pathogen population
- **C**: Commensal population
- **E**: Epithelial health score (0 = fully damaged, E_max = fully healthy)
- **ROS**: Reactive oxygen species concentration
- **DAMPS**: Damage-associated molecular patterns
- **PAMPS**: Pathogen-associated molecular patterns
- **SAMPS**: Suppression-associated molecular patterns

## Key Files

### Core ODE System
- **`MISC/ODE_SYSTEM.R`** - Main ODE implementation
  - `ode_system()`: Differential equation system
  - `setup_initial_conditions()`: Initialize state variables
  - `run_ode_simulation()`: ODE solver wrapper
  - `process_ode_output()`: Convert ODE output to ABM-compatible format

### Data Generation Scripts
- **`DLL_datagen_ode.R`** - Main data generation orchestrator (analogous to `DLL_datagen_abm.R`)
  - Reads parameter sets from `lhs_parameters_della.csv`
  - Loops over scenarios and parameter sets
  - Saves results to `/scratch/gpfs/CMETCALF/sim_ode/`

- **`MISC/RUN_REPS_ODE_PAMPS.R`** - Run ODE for multiple replicates (analogous to `RUN_REPS_CPP_ABM_PAMPS.R`)
  - Sets up ODE parameters
  - Runs ODE solver
  - Processes output

### Testing & Documentation
- **`test_ode_single_param.R`** - Test script for single parameter set
  - Validates ODE implementation
  - Generates diagnostic plots
  - Outputs summary statistics

- **`ODE_PARAMETER_MAPPING.md`** - Detailed parameter mapping guide
  - ABM → ODE parameter translations
  - New ODE-specific parameters
  - Expected differences between ABM and ODE

## Usage

### Quick Test (Single Parameter Set)
```bash
Rscript test_ode_single_param.R
```

This runs the ODE system with parameter set 0, pathogen level 2, for 3 replicates, and generates:
- `test_ode_output.rds` - Full longitudinal data
- `test_ode_populations.png` - Population dynamics plot
- `test_ode_epithelial.png` - Epithelial health plot
- `test_ode_signals.png` - Signaling molecule dynamics

### Full Data Generation (All Parameter Sets)
```bash
# Run on cluster with parallel jobs
# Split into N chunks, run chunk M:
Rscript DLL_datagen_ode.R <N> <M>

# Example: Split into 100 jobs, run job 1:
Rscript DLL_datagen_ode.R 100 1
```

Output files saved to `/scratch/gpfs/CMETCALF/sim_ode/`:
```
longitudinal_df_param_set_id_<ID>_sterile_<0/1>_macspec_<0/1/2>_tregs_<0/1>_ros_level_<LEVEL>_pat_level_<LEVEL>_trnd_<0/1>.rds
```

## Key Mechanisms

### 1. Macrophage Polarization
```
M0 → M1 when: danger_signal > safety_signal AND danger_signal > threshold
M0 → M2 when: safety_signal > danger_signal AND safety_signal > threshold

danger_signal = DAMPS + PAMPS
safety_signal = SAMPS
```

Polarization rate uses Michaelis-Menten kinetics:
```
rate = polarization_speed × M0 × (signal / (signal + K_polarization))
```

### 2. Epithelial Dynamics
```
dE/dt = recovery - damage_from_pathogens - damage_from_ROS

recovery = epith_recovery_chance × (E_max - E)
damage_from_pathogens = damage_rate_pathogen × P × E / (K_pathogen_damage + P)
damage_from_ROS = damage_rate_ROS × (ROS - threshold) × E / E_max   [if ROS > threshold]
```

### 3. Microbe Dynamics
```
dP/dt = leakage - ROS_killing - phagocytosis

leakage = rate_leak_pathogen_injury × injury_fraction × grid_size
ROS_killing = kill_rate_ROS × (ROS - threshold) × P   [if ROS > threshold]
phagocytosis = total_phagocytic_capacity × (P / (P + C)) × P / (K_phago + P + C)
```

### 4. Signaling Molecules
```
dROS/dt = activity_ROS_M1 × M1 × add_ROS - ros_decay × ROS
dDAMPS/dt = injury_production + microbe_production - DAMPs_decay × DAMPS
dPAMPS/dt = add_PAMPs × P - PAMPs_decay × PAMPS
dSAMPS/dt = treg_activation × add_SAMPs - SAMPs_decay × SAMPS
```

## ODE-Specific Parameters

These parameters are unique to the ODE version and need calibration:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `polarization_speed` | 0.1 | Rate of M0→M1/M2 conversion |
| `K_polarization` | 0.5 | Half-saturation for polarization |
| `damage_rate_pathogen` | 1.0 | Pathogen-induced epithelial damage rate |
| `damage_rate_ROS` | 2.0 | ROS-induced epithelial damage rate |
| `K_pathogen_damage` | 50 | Half-saturation for pathogen damage |
| `kill_rate_ROS` | 2.0 | ROS killing efficiency |
| `K_phagocytosis` | 100 | Half-saturation for phagocytosis |
| `treg_activation_efficiency` | 0.01 | Treg activation rate |
| `K_treg_activation` | 50 | Half-saturation for Treg activation |
| `K_DAMPS_microbe` | 100 | Half-saturation for microbe-induced DAMPs |

## Output Format

The ODE system produces output in the same format as the ABM for easy comparison:

### Columns
- `t` - Time step (0 to t_max)
- `M0, M1, M2` → `phagocyte_M0, phagocyte_M1, phagocyte_M2`
- `P` → `pathogen`
- `C` → `commensal`
- `E` → `epithelial_score`
- `ROS, DAMPS, PAMPS, SAMPS` - Signal concentrations
- `C_ROS, C_M1, C_M2` - Cumulative commensal deaths (by cause)
- `P_ROS, P_M1, P_M2` - Cumulative pathogen deaths (by cause)
- `time_ss_e, time_ss_p` - Time to steady state (epithelial, pathogen)
- Scenario metadata: `param_set_id, rep_id, sterile, macspec_on, tregs_on, ros_level, pat_level, randomize_tregs`

### Note on Replicates
Since the ODE is deterministic, all replicates with the same parameters produce identical results. We keep `num_reps` for consistency with ABM output structure, but ODE "replicates" are exact duplicates.

## Comparison with ABM

### Advantages of ODE Version
✓ **Much faster**: No spatial loops, no stochastic events
✓ **Deterministic**: Reproducible results
✓ **Smooth dynamics**: No discrete event noise
✓ **Easier analysis**: Standard ODE theory applies

### Limitations vs ABM
✗ **No spatial structure**: Can't capture local effects (ROS hotspots, clustering)
✗ **No stochasticity**: No variability between replicates
✗ **No agent memory**: Can't implement bacteria registry or discrimination
✗ **No discrete events**: Individual engulfment, injury events lost
✗ **Mean-field approximation**: Assumes well-mixed, ignores correlations

### Expected Differences
The ODE version should:
1. Show similar **qualitative behavior** (e.g., pathogen clearance or persistence)
2. Have similar **steady-state values** (with calibration)
3. Differ in **transient dynamics** (smoother, no fluctuations)
4. Miss **spatial phenomena** (e.g., local Treg-macrophage interactions)

## Calibration Strategy

To match ABM behavior:

1. **Select representative parameter sets** (e.g., 10-20 spanning parameter space)
2. **Run both ABM and ODE** with same parameters
3. **Compare key metrics**:
   - Final epithelial score
   - Final pathogen count
   - Time to steady state
   - Peak pathogen load
   - M1:M2 ratio over time
4. **Adjust ODE-specific parameters** to minimize differences
5. **Validate** on held-out parameter sets

### Suggested Metrics for Calibration
```r
# Root mean squared error for epithelial score trajectory
rmse_epithelial = sqrt(mean((ode$epithelial_score - abm$epithelial_score)^2))

# Final state error
error_final_pathogen = abs(ode$pathogen[t_max] - mean(abm$pathogen[t_max]))

# Steady state time difference
error_ss_time = abs(ode$time_ss_e - mean(abm$time_ss_e))
```

## Solver Details

The ODE system uses R's `deSolve` package with the `lsoda` solver:
- **Adaptive step size**: Automatically adjusts for stiff equations
- **Error control**: Maintains accuracy while optimizing speed
- **Stiffness handling**: Switches between methods for stiff/non-stiff regions

Typical solve time: ~1-2 seconds per simulation (vs ~10-60 seconds for ABM)

## Dependencies

Required R packages:
```r
install.packages(c("deSolve", "dplyr", "tidyr", "zoo", "ggplot2"))
```

## Future Extensions

Possible enhancements:
1. **Stochastic differential equations (SDEs)**: Add noise to capture variability
2. **Spatial structure**: Partial differential equations (PDEs) with diffusion
3. **Age structure**: Macrophage age classes for memory effects
4. **Parameter inference**: Fit ODE parameters to ABM data
5. **Sensitivity analysis**: Fast exploration of parameter space
6. **Bifurcation analysis**: Identify critical parameter thresholds

## Contact & Citation

This ODE version is a mean-field approximation of the spatial ABM described in:
- Original ABM documentation: `Agent-Based Model Simulation Rules.md`
- Parameter generation: `DLL_paramgen_abm.R`
- Full ABM implementation: `MISC/RUN_REPS_CPP_ABM_PAMPS.R`

For questions or issues with the ODE implementation, see the parameter mapping guide (`ODE_PARAMETER_MAPPING.md`) or run the test script (`test_ode_single_param.R`) for diagnostics.
