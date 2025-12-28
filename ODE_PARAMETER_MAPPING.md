# ODE Parameter Mapping

This document explains how parameters from the agent-based model (ABM) are mapped to the ODE system.

## Direct Mappings (Unchanged)

These parameters are used directly in the ODE system with the same interpretation:

### Thresholds
- `th_ROS_microbe` - ROS concentration threshold for microbe killing
- `th_ROS_epith_injury` - ROS concentration threshold for epithelial injury
- `epith_recovery_chance` - Epithelial recovery rate (per timestep)

### Activation Thresholds
- `activation_threshold_danger` - Danger signal threshold for M1 polarization
- `activation_threshold_SAMPs` - SAMPs threshold for M2 polarization

### Macrophage Activities
- `activity_ROS_M0_baseline`, `activity_ROS_M1_baseline`, `activity_ROS_M2_baseline` - ROS production rates
- `activity_engulf_M0_baseline`, `activity_engulf_M1_baseline`, `activity_engulf_M2_baseline` - Phagocytosis rates

### Microbe Leakage Rates
- `rate_leak_pathogen_injury` - Pathogen leakage rate (injury-dependent)
- `rate_leak_commensal_baseline` - Baseline commensal leakage
- `rate_leak_commensal_injury` - Injury-enhanced commensal leakage

### Recruitment
- `recruitment_rate_danger` - Danger-dependent macrophage recruitment

### Signaling Molecule Parameters
- `add_ROS`, `add_DAMPs`, `add_PAMPs`, `add_SAMPs` - Production rates
- `ros_decay`, `DAMPs_decay`, `PAMPs_decay`, `SAMPs_decay` - Decay rates

### Environment Parameters
- `grid_size` - Used to scale population-level dynamics
- `n_phagocytes` - Initial macrophage population (all M0)
- `n_commensals_lp` - Initial commensal population
- `injury_percentage` - Initial epithelial damage percentage
- `max_level_injury` - Maximum injury level (5 in ABM)

## Parameters NOT Used in ODE

These ABM-specific parameters are not needed in the well-mixed ODE system:

### Spatial Parameters (No Spatial Structure in ODE)
- `diffusion_speed_*` - All diffusion speeds (no spatial structure to diffuse through)
- `act_radius_*` - Activation radii (no local neighborhoods)

### Agent-Specific Parameters
- `cc_phagocyte` - Bacteria registry capacity (ODE doesn't track individual macrophage memory)
- `active_age_limit` - Phenotype reassessment time (ODE uses continuous polarization)
- `treg_vicinity_effect` - Treg activation radius (no spatial structure)
- `treg_discrimination_efficiency` - Commensal/pathogen discrimination (simplified in ODE)
- `mac_discrimination_efficiency` - Macrophage discrimination (simplified in ODE)
- `macspec_on` - Macrophage specificity mode (not implemented in basic ODE version)
- `randomize_tregs` - Treg movement mode (no movement in ODE)

### Stochastic/Movement Parameters
- All parameters related to agent movement and spatial positioning
- Stochastic recognition parameters

## New ODE-Specific Parameters

These parameters are introduced for the ODE system and need to be calibrated:

### Polarization Dynamics
- `polarization_speed = 0.1` - Rate constant for M0→M1/M2 conversion
  - Higher values = faster polarization response
  - Default 0.1 means ~10% of M0s polarize per timestep when signals are saturating
- `K_polarization = 0.5` - Half-saturation constant for polarization
  - Signal level at which polarization rate is half-maximal
  - Represents threshold between slow and fast polarization

### Epithelial Damage
- `damage_rate_pathogen = 1.0` - Rate of pathogen-induced epithelial damage
  - Per pathogen per timestep, saturating with Michaelis-Menten kinetics
- `damage_rate_ROS = 2.0` - Rate of ROS-induced epithelial damage
  - Linear with ROS above threshold
- `K_pathogen_damage = 50` - Half-saturation constant for pathogen damage
  - Number of pathogens at which damage rate is half-maximal

### Microbe Killing
- `kill_rate_ROS = 2.0` - ROS killing efficiency (per unit ROS above threshold)
  - Higher values = more effective ROS killing
- `K_phagocytosis = 100` - Half-saturation constant for phagocytosis
  - Total microbe count at which phagocytosis rate is half-maximal
  - Represents saturation of phagocytic capacity

### Treg Dynamics
- `treg_activation_efficiency = 0.01` - Rate of Treg activation
  - Per macrophage per commensal per timestep
- `K_treg_activation = 50` - Half-saturation for Treg activation
  - Commensal count at which activation rate is half-maximal

### Signaling
- `K_DAMPS_microbe = 100` - Half-saturation for microbe-induced DAMPs
  - Microbe count at which DAMP production is half-maximal

## Epithelial Health Score

In the ODE system:
- `E` represents total epithelial health score (continuous variable)
- `E_max = grid_size^2 * 6` = fully healthy epithelium (all 625 cells at health level 6)
- `E = 0` = completely damaged epithelium
- Injury fraction = `1 - (E / E_max)`, ranges from 0 (healthy) to 1 (fully damaged)

In the ABM:
- Epithelium tracked as discrete cells with injury levels 0-5
- Health score = sum over cells: `6×(level 0) + 5×(level 1) + ... + 1×(level 5)`

The ODE version simplifies this to a single continuous health score.

## Compartment Correspondences

| ODE Compartment | ABM Equivalent | Notes |
|-----------------|----------------|-------|
| `M0, M1, M2` | `phagocyte_phenotype` counts | Direct correspondence |
| `P` | `pathogen_coords` rows | Population counts |
| `C` | `commensal_coords` rows | Population counts |
| `E` | `epithelial_score` | Continuous vs discrete |
| `ROS, DAMPS, PAMPS, SAMPS` | Signal matrices (summed) | Well-mixed vs spatial |

## Calibration Strategy

To match ABM behavior, the ODE-specific parameters should be calibrated by:

1. Run ABM with a parameter set
2. Run ODE with same parameter set (using default ODE-specific values)
3. Compare key outputs:
   - Epithelial health trajectories
   - Pathogen/commensal dynamics
   - Macrophage polarization ratios
   - Time to steady state
4. Adjust ODE-specific parameters to minimize differences
5. Repeat across multiple parameter sets

Key metrics for comparison:
- Final steady-state values (epithelial score, pathogen count)
- Time to steady state (`time_ss_e`, `time_ss_p`)
- Peak pathogen load
- M1:M2 ratio dynamics

## Expected Differences

Even with perfect calibration, some differences are inherent:

1. **No stochasticity**: ODE is deterministic
   - Solution: Run multiple "replicates" with slight parameter variations
2. **No spatial heterogeneity**: Can't capture local effects
   - Local ROS hotspots killing nearby microbes
   - Spatial clustering of macrophages
   - Treg-macrophage co-localization
3. **No discrete events**: Everything is continuous
   - Individual engulfment events become continuous rates
   - Individual cell injury becomes continuous damage
4. **Simplified discrimination**: No agent-level memory
   - Can't implement bacteria registry or two-step recognition

The ODE version is best viewed as a **mean-field approximation** of the ABM, capturing population-level dynamics while losing individual-level and spatial details.
