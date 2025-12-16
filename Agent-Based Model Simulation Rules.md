# Agent-Based Model Simulation Rules (Updated)

## Overview
This document describes the agent-based model (ABM) implemented in `DLL_datagen_abm.R` and associated scripts in the `./MISC` directory. The simulation models immune responses in the gut lamina propria following epithelial injury, including interactions between epithelial cells, microbes (commensals and pathogens), phagocytes (macrophages/dendritic cells), and regulatory T cells (Tregs).

**Key Features:**
- Spatial simulation on a 2D grid
- Multiple signaling molecules (DAMPs, PAMPs, SAMPs, ROS)
- Phagocyte polarization (M0, M1, M2)
- Treg activation and suppression
- Macrophage specificity (discrimination between commensals and pathogens)
- C++ acceleration for computational efficiency (20-100x speedup)

---

## 1. General Parameters & Environment

### Grid and Time
- **Grid Size:** The simulation takes place on a 25×25 grid (`grid_size = 25`)
- **Time Steps:** The simulation runs for a maximum of 5000 time steps (`t_max = 5000`)
- **Replicates:** Each parameter set is simulated 10 times (`num_reps = 10`) to account for stochasticity

### Initial Populations
- **Phagocytes:** `n_phagocytes = round(grid_size * grid_size * 0.20)` (125 cells)
- **Tregs:** `n_tregs = round(grid_size * grid_size * 0.20)` (125 cells)
- **Commensals (initial):** `n_commensals_lp = 20`
- **Pathogens (initial):** Depends on injury site size and leak rate

### Injury Definition
- **Injury Site:** Defined as the middle 60% of the epithelial layer (`injury_percentage = 60`)
- **Max Injury Level:** `max_level_injury = 5`
- **Injury Dynamics:** Determined by logistic function parameters `k_in = 0.044` and `x0_in = 50`

### Signal Constraints
All signaling molecules have maximum values to prevent unbounded accumulation:
- `max_cell_value_ROS = 1`
- `max_cell_value_DAMPs = 1`
- `max_cell_value_SAMPs = 1`
- `max_cell_value_PAMPs = 1`

### Activation Radii
- `act_radius_ROS = 1`
- `act_radius_treg = 1`
- `act_radius_DAMPs = 1`
- `act_radius_SAMPs = 1`

### Stochasticity & Classification
Since the same parameter set can yield different outcomes depending on the random seed, each set is simulated 10 times. A parameter set is assigned to a class—Tregs helping, harming, or being irrelevant to epithelial health or pathogen control—if the corresponding outcome is observed more than half the time across the realizations (more than 5 runs in 10 runs per parameter set).

---

## 2. Epithelial Cells

### Location & Structure
- **Location:** Epithelial cells are located at y=0 (effectively y=1 in the 1-indexed grid for interaction purposes)
- **Grid Representation:** Each x-coordinate (1 to `grid_size`) represents one epithelial cell

### Injury Levels
Epithelial cells have 6 states:
- **0:** Healthy
- **1-5:** Injured, with 5 being the maximum injury level

### Initial Condition
The simulation begins with epithelial cells within the `injury_site` at injury level 1.

### Injury Progression
Epithelial cells can become more injured through three mechanisms:

1. **Pathogen-induced injury:**
   - Pathogens at the epithelium (y=1) directly increase injury levels
   - Injury increment calculated using `logistic_scaled_0_to_5_quantized(pathogen_count)`
   - Based on logistic function that maps pathogen count to injury level (0-5)

2. **ROS-induced injury:**
   - If mean ROS level in the vicinity (radius `act_radius_ROS`) exceeds threshold `th_ROS_epith_injury`, injury level increases by 1
   - Calculated using C++ accelerated function `calculate_epithelial_ros_cpp()`

3. **Maximum constraint:**
   - Injury level cannot exceed `max_level_injury = 5`

### Recovery
- **Stochastic recovery:** At each time step, injured epithelial cells have a probability `epith_recovery_chance` of reducing their injury level by 1
- Recovery is processed cell-by-cell to maintain consistent random number consumption

### DAMP Release
Epithelial cells release DAMPs (Damage-Associated Molecular Patterns) through two mechanisms:

1. **Injury-proportional release:**
   - DAMPs released proportional to `level_injury * add_DAMPs`
   - Applied at y=1 in the grid

2. **Microbe detection at epithelium:**
   - Both commensals and pathogens at the epithelium (y=1) trigger DAMP release
   - The epithelium cannot differentiate between commensals and pathogens based on this PRR (Pattern Recognition Receptor) stimulation
   - Additional DAMPs calculated using `logistic_scaled_0_to_5_quantized(pathogen_count + commensal_count) * add_DAMPs`
   - Represents basolateral PRR stimulation by microbes

### Microbe Leakage
The epithelial barrier integrity decreases with injury:

- **Pathogen leakage:**
  - Rate: `rate_leak_pathogen_injury * mean(level_injury) * length(injury_site_updated)`
  - New pathogens leak at each time step proportional to epithelial injury
  - If `sterile = 1`, no new pathogens leak (sterile injury scenario)
  - Pathogens preferentially leak through more injured sites (sampling weighted by `level_injury`)

- **Commensal leakage:**
  - Baseline rate: `rate_leak_commensal_baseline * grid_size` (across entire epithelium)
  - Injury-enhanced rate: `rate_leak_commensal_injury * mean(level_injury) * length(injury_site_updated)`
  - Total commensals entering: baseline + injury-enhanced
  - Injury-enhanced commensals preferentially leak through injured sites

---

## 3. Microbes

### Types
- **Commensals:** Non-pathogenic bacteria that normally colonize the gut
- **Pathogens:** Virulent bacteria that damage epithelial cells and produce PAMPs

### Movement
Microbes perform biased random walks:

1. **At epithelium (y=1):**
   - Can only move deeper into tissue (dy ∈ {1}) or laterally/stay (if selected)
   - Horizontal movement (dx) is random when dy=0
   - This prevents microbes from leaving the lamina propria once they enter

2. **Within lamina propria (y>1):**
   - Can move in any of 8 directions or stay still
   - Movement is uncorrelated and random
   - Bounded by grid limits (1 to `grid_size`)

3. **Implementation:**
   - **IMPORTANT:** Microbes move BEFORE DAMP calculation (updated order)
   - This ensures DAMPs reflect current microbe positions

### Leakage into Lamina Propria
(See Epithelial Cells - Microbe Leakage section above)

### Death Mechanisms

1. **ROS-induced death:**
   - Microbes die if local ROS concentration exceeds `th_ROS_microbe`
   - Death check performed over a radius `act_radius_ROS` around each microbe
   - Implemented using C++ accelerated function `kill_microbes_with_ros_cpp()` for 50-100x speedup
   - Separate death counters tracked for commensals and pathogens

2. **Phagocytosis:**
   - Microbes are engulfed by phagocytes at the same location
   - Probability of engulfment depends on phagocyte's `activity_engulf`
   - Tracked separately by phagocyte phenotype (M0, M1, M2)

### Signal Production

**Pathogens produce PAMPs (Pathogen-Associated Molecular Patterns):**
- PAMPs are toxins released by pathogens that diffuse through the environment
- Production rate: `add_PAMPs * logistic_scaled_0_to_5_quantized(pathogen_count)` at each pathogen location
- PAMPs provide an additional danger signal distinct from DAMPs
- **Key difference from DAMPs:** PAMPs are specific to pathogens, while DAMPs indicate general tissue damage

---

## 4. Phagocytes (Macrophages/Dendritic Cells)

### Population
- Initial number: `n_phagocytes = round(grid_size * grid_size * 0.20)` = 125 cells
- Do not proliferate or die during simulation
- **NEW:** Can be dynamically recruited from tissue borders in response to danger signals (see Recruitment section below)

### Movement & Chemotaxis
Phagocytes perform chemotaxis toward danger signals:

- **Chemotaxis gradient:** `density_matrix_phagocytes = DAMPs + PAMPs`
  - **KEY UPDATE:** Phagocytes respond to combined DAMPs and PAMPs
  - This allows them to respond to both tissue damage (DAMPs) and pathogen signals (PAMPs)

- **Movement algorithm:**
  - If gradient exists: Move probabilistically toward higher combined signal (DAMPs + PAMPs)
  - If no gradient (all values equal): Move randomly (biased away from y=0 boundary)
  - Uses 3×3 neighborhood (Moore neighborhood)

### Phenotypes
Phagocytes can exist in three phenotypes:

1. **M0 (Resting):** Default state
   - Minimal ROS production: `activity_ROS_M0_baseline = 0`
   - Low engulfment activity: `activity_engulf_M0_baseline = 0.05`

2. **M1 (Pro-inflammatory):** Activated by danger signals
   - High ROS production: `activity_ROS_M1_baseline` (from parameter set)
   - High engulfment activity: `activity_engulf_M1_baseline` (from parameter set)

3. **M2 (Anti-inflammatory/Resolving):** Activated by safety signals
   - No ROS production: `activity_ROS_M2_baseline = 0`
   - Moderate engulfment activity: `activity_engulf_M2_baseline` (from parameter set)

### Activation Logic (M0 → M1 or M2)

M0 phagocytes assess their local environment and differentiate based on signal dominance:

**Signal Calculation:**
- `danger_signal = avg_DAMPs + avg_PAMPs` (combined danger signal)
- `avg_SAMPs` calculated over radius `act_radius_SAMPs`
- Calculated using C++ accelerated function `calculate_phagocyte_signals_cpp()`

**Differentiation Rules:**

1. **M0 → M1 (Pro-inflammatory):**
   - Condition: `danger_signal >= activation_threshold_danger AND danger_signal > avg_SAMPs`
   - Initiates when danger signals dominate over safety signals
   - Sets `phagocyte_active_age = 1`

2. **M0 → M2 (Anti-inflammatory):**
   - Condition: `avg_SAMPs >= activation_threshold_SAMPs AND avg_SAMPs > danger_signal`
   - Initiates when safety signals dominate over danger signals
   - Sets `phagocyte_active_age = 1`

### Phenotype Plasticity (M1 ↔ M2 ↔ M0)

After being active for `active_age_limit` time steps, M1 and M2 phagocytes reassess their environment:

**Without Macrophage Specificity (`macspec_on = 0`):**
- Re-evaluates signals as during initial activation
- Can switch: M1 → M2, M2 → M1, or either → M0
- Based purely on current environmental signals

**With Macrophage Specificity (`macspec_on > 0`):**
Phagocytes use both environmental signals AND engulfment history:

1. **Calculate engulfment pattern:**
   - `num_pat_engulfed` = count of pathogens in bacteria registry
   - `num_com_engulfed` = count of commensals in bacteria registry
   - `rat_com_pat_real` = true ratio of commensal to total bacteria

2. **Apply discrimination:**
   - **If `macspec_on = 1`:** Use same discrimination as Tregs
     - `mac_discrimination_efficiency = treg_discrimination_efficiency`
   - **If `macspec_on = 2`:** Perfect discrimination
     - `mac_discrimination_efficiency = 1.0`

   - Perceived ratio: `rat_com_pat = mac_discrimination_efficiency * rat_com_pat_real + (1 - mac_discrimination_efficiency) * runif(1)`
   - Interpolates between truth (discrimination=1) and random guess (discrimination=0)

3. **Determine dominance:**
   - Environmental: `DAMPs_dominant` vs `SAMPs_dominant` (using danger = DAMPs + PAMPs)
   - Engulfment: `pathogen_engulfment_dominant` (rat_com_pat ≤ 1 - threshold) vs `commensal_engulfment_dominant` (rat_com_pat > threshold)

4. **Polarization logic (danger-biased):**
   - **→ M1:** If `DAMPs_dominant OR pathogen_engulfment_dominant`
     - Either environmental danger OR pathogen engulfment is sufficient
   - **→ M2:** If `SAMPs_dominant AND commensal_engulfment_dominant`
     - BOTH environmental safety AND commensal engulfment required
   - **→ M0:** If both signals below threshold

### Bacteria Registry & Digestion
Each phagocyte maintains a memory of recently engulfed bacteria:

- **Registry:** `phagocyte_bacteria_registry` matrix (size `n_phagocytes × cc_phagocyte`)
  - Entries: `+1` = commensal, `-1` = pathogen, `0` = empty
  - Organized as a shift register (newest entries at the beginning)

- **Capacity:** `cc_phagocyte` = carrying capacity (max bacteria processed simultaneously)
  - From parameter set (typically 20-30)

- **Insertion:** New engulfments added to front using C++ function `shift_insert_fast_cpp()`
  - Oldest entries pushed out when capacity reached

- **Usage:** Used for calculating discrimination-based polarization (when `macspec_on > 0`)

### ROS Production
- Only M1 phagocytes produce ROS
- ROS added at phagocyte location: `ROS[y, x] += activity_ROS * add_ROS`
- No direct radius effect (ROS diffuses after production)

### Engulfment Process
At each time step, for each phagocyte:

1. **Find collocated microbes:**
   - Check pathogens and commensals at same (x, y) location

2. **Attempt engulfment:**
   - For each microbe: Success if `runif(1) < activity_engulf`
   - Success rate depends on phenotype (M1 and M2 higher than M0)

3. **Update registry:**
   - Successful engulfments added to `bacteria_registry`
   - Microbe removed from environment

4. **Track deaths:**
   - Separate counters for commensals/pathogens killed by M0/M1/M2

### Recruitment from Tissue Borders

**NEW FEATURE - Border Recruitment Based on Local Danger Signals**

Macrophages can be dynamically recruited into the tissue from the three non-epithelial borders in response to inflammation:

**Recruitment Mechanism:**
1. **Recruitment borders:**
   - **TOP border** (y = grid_size): All x positions
   - **LEFT border** (x = 1): All y positions except epithelium (y > 1)
   - **RIGHT border** (x = grid_size): All y positions except epithelium (y > 1)

2. **Location-specific recruitment:**
   - For each border position, calculate local danger signal: `danger_at_border = DAMPs[border_pos] + PAMPs[border_pos]`
   - Number recruited: `n_recruit ~ Poisson(λ = recruitment_rate_danger × danger_at_border)`
   - Stochastic recruitment proportional to local inflammation

3. **Properties of recruited macrophages:**
   - Enter at exact border position with high danger signal
   - Initialize as M0 phenotype (resting/naive)
   - Have empty bacteria registry
   - Immediately begin sensing local signals and can activate

4. **Timing:**
   - Recruitment occurs after signal diffusion/decay
   - Before cell movement
   - Ensures recruited cells respond to current danger landscape

**Parameters:**
- `recruitment_rate_danger` (range: 0-0.5)
  - Controls recruitment rate via Poisson λ parameter
  - Higher values = more aggressive recruitment
  - Set to 0 disables recruitment (backward compatibility)

**Biological Interpretation:**
- Mimics monocyte recruitment from blood vessels at tissue periphery
- Recruitment strength proportional to local chemokine gradients
- Realistic spatial heterogeneity in immune cell infiltration
- Allows model to capture dynamic changes in immune cell numbers during inflammation

---

## 5. Tregs (Regulatory T Cells)

### Population
- Fixed number: `n_tregs = round(grid_size * grid_size * 0.20)` = 125 cells
- Do not proliferate or die during simulation

### Movement & Chemotaxis

Treg movement depends on the `randomize_tregs` scenario parameter:

- **Normal mode (`randomize_tregs = 0`):**
  - `density_matrix_tregs = DAMPs`
  - Tregs chemotax toward DAMPs (sites of tissue damage)
  - Uses same chemotaxis algorithm as phagocytes

- **Random mode (`randomize_tregs = 1`):**
  - `density_matrix_tregs = 0 * DAMPs` (all zeros)
  - Tregs move randomly (no chemotaxis)
  - Used as a control to test importance of Treg localization

### Phenotypes
Tregs have two states:

1. **Resting (phenotype = 0):**
   - Default state
   - No SAMP production
   - `activity_SAMPs_binary = 0`

2. **Activated (phenotype = 1):**
   - Produces SAMPs
   - `activity_SAMPs_binary = 1`

### Activation Logic

Tregs are activated by encountering phagocytes presenting predominantly commensal antigens:

**Prerequisites:**
- `allow_tregs = 1` (Tregs functional in this scenario)
- Treg within `treg_vicinity_effect` distance of M1 or M2 phagocyte
- Phagocyte has engulfed at least one bacterium

**Antigen Recognition Process:**

1. **Calculate antigen ratio:**
   - `num_pat_antigens` = pathogens in phagocyte bacteria registry
   - `num_com_antigens` = commensals in phagocyte bacteria registry
   - `rat_com_pat_real` = true commensal ratio

2. **Apply discrimination:**
   - Tregs have imperfect discrimination ability: `treg_discrimination_efficiency`
   - Perceived ratio: `rat_com_pat = treg_discrimination_efficiency * rat_com_pat_real + (1 - treg_discrimination_efficiency) * runif(1)`
   - When discrimination = 0: completely random guess
   - When discrimination = 1: perfect recognition
   - Values between 0-1: interpolation between random and perfect

3. **Activation decision:**
   - If `rat_com_pat > rat_com_pat_threshold`: Treg activates
   - Threshold typically 0.5-0.9 (from parameter set)
   - All nearby Tregs (within `treg_vicinity_effect`) activate simultaneously

**Upon activation:**
- `treg_phenotype = 1`
- `treg_activity_SAMPs_binary = 1`
- `treg_active_age = 1`

### SAMP Production
- Activated Tregs produce SAMPs at their location
- Production: `SAMPs[y, x] += add_SAMPs * treg_activity_SAMPs_binary`
- Uses C++ accelerated function `update_SAMPs_batch_cpp()` for all active Tregs

### Deactivation
- After `active_age_limit` time steps in activated state, Tregs return to resting
- All parameters reset to 0

---

## 6. Signaling Molecules (DAMPs, PAMPs, SAMPs, ROS)

### DAMPs (Damage-Associated Molecular Patterns)

**Sources:**
1. Injured epithelial cells (proportional to injury level)
2. Microbes at epithelium triggering PRR stimulation
3. (Legacy code: pathogens within lamina propria, but multiplied by 0)

**Parameters:**
- Production rate: `add_DAMPs`
- Diffusion speed: `diffusion_speed_DAMPs`
- Decay rate: `DAMPs_decay`
- Maximum value: `max_cell_value_DAMPs = 1`

**Function:**
- Signals tissue damage and danger
- Activates phagocytes → M1
- Attracts phagocytes and (optionally) Tregs

### PAMPs (Pathogen-Associated Molecular Patterns)

**NEW FEATURE - Not in original documentation**

**Sources:**
- Released by pathogens at their location
- Production proportional to pathogen count: `add_PAMPs * logistic_scaled_0_to_5_quantized(pathogen_count)`

**Parameters:**
- Production rate: `add_PAMPs`
- Diffusion speed: `diffusion_speed_PAMPs`
- Decay rate: `PAMPs_decay`
- Maximum value: `max_cell_value_PAMPs = 1`

**Function:**
- Signals pathogen presence (toxins)
- Combines with DAMPs to form "danger signal" for phagocyte activation
- Attracts phagocytes (chemotaxis toward DAMPs + PAMPs)
- Key innovation: Provides pathogen-specific signal distinct from tissue damage

**Biological interpretation:**
- Represents pathogen toxins (e.g., LPS, flagellin, peptidoglycans)
- Diffuses through tissue similar to DAMPs
- Allows immune system to distinguish pathogenic infection from sterile injury

### SAMPs (Suppression-Associated Molecular Patterns)

**Sources:**
- Released by activated Tregs

**Parameters:**
- Production rate: `add_SAMPs`
- Diffusion speed: `diffusion_speed_SAMPs`
- Decay rate: `SAMPs_decay`
- Maximum value: `max_cell_value_SAMPs = 1`

**Function:**
- Signals immunosuppression and tissue repair
- Activates phagocytes → M2
- Counter-balances danger signals (DAMPs + PAMPs)

**Biological interpretation:**
- Represents anti-inflammatory cytokines (e.g., IL-10, TGF-β)

### ROS (Reactive Oxygen Species)

**Sources:**
- Released by M1 phagocytes
- Production: `activity_ROS_M1_baseline * add_ROS`

**Parameters:**
- Production rate: `add_ROS`
- Diffusion speed: `diffusion_speed_ROS`
- Decay rate: `ros_decay`
- Maximum value: `max_cell_value_ROS = 1`

**Function:**
- Antimicrobial activity (kills microbes)
- Causes epithelial injury (collateral damage)

**Control scenario:**
- When `control = 1`: `add_ROS = 0` and `activity_ROS_M1_baseline = 0`
- Tests whether infection can resolve without ROS

### Diffusion & Decay

**Diffusion:**
- All signals diffuse using discrete Laplacian operator
- Implemented with C++ function `diffuse_matrix_cpp()` for 5-10x speedup
- Boundary conditions: Reflective on all sides except top (y=0 has no reflection)
- Diffusion speeds are parameters sampled from parameter set

**Decay:**
- Linear decay at each time step: `Signal = Signal - decay_rate * Signal`
- Prevents unbounded accumulation
- Different decay rates for each signal type

**Maximum values:**
- All signals capped at maximum value (typically 1.0)
- Prevents numerical instability

---

## 7. Simulation Algorithm (Order of Operations)

At each time step `t`, the following operations are performed in this exact order:

### 1. Update Injury Site
- Identify currently injured epithelial cells: `injury_site_updated = which(level_injury > 0)`

### 2. Update SAMPs (from activated Tregs)
- For each activated Treg (phenotype = 1):
  - Add SAMPs at Treg location
  - Uses C++ batch function `update_SAMPs_batch_cpp()`

### 3. Update ROS (from M1 phagocytes)
- For each M1 phagocyte:
  - Add ROS at phagocyte location: `ROS[y, x] += activity_ROS * add_ROS`

### 4. Move Microbes
- **Pathogens and commensals** perform random walks:
  - At epithelium (y=1): Can only move up or laterally
  - In lamina propria (y>1): Can move in any direction
  - Coordinates clamped to grid boundaries [1, grid_size]

### 5. Pre-Calculate Microbe Counts at Epithelium
- Count pathogens and commensals at y=1 for each x position
- Used for DAMP generation and epithelial injury calculation

### 6. Update DAMPs
- Add DAMPs from injured epithelium: `DAMPs[1, x] += level_injury[x] * add_DAMPs`
- Add DAMPs from microbes at epithelium: `DAMPs[1, x] += logistic_scaled(pathogen_count + commensal_count) * add_DAMPs`

### 7. Update PAMPs
- For each grid position with pathogens:
  - Add PAMPs: `PAMPs[y, x] += logistic_scaled(pathogen_count) * add_PAMPs`

### 8. Diffuse & Decay Signals
- **Diffusion:** All signals (DAMPs, PAMPs, SAMPs, ROS) diffuse via discrete Laplacian
  - Uses C++ function `diffuse_matrix_cpp()` for 5-10x speedup
  - Reflective boundary conditions (except top)
- **Decay:** Linear decay for all signals: `Signal = Signal * (1 - decay_rate)`
- **Capping:** Enforce maximum values for all signals

### 9. Recruit Macrophages from Borders (NEW)
- If `recruitment_rate_danger > 0`:
  - For each border position (top, left, right):
    - Calculate local danger: `danger_at_border = DAMPs + PAMPs`
    - Recruit: `n_recruit ~ Poisson(recruitment_rate_danger * danger_at_border)`
    - Add new M0 macrophages at border position

### 10. Move Phagocytes and Tregs
- **Phagocytes:** Chemotax toward `DAMPs + PAMPs`
- **Tregs:** Chemotax toward `DAMPs` (if `randomize_tregs = 0`) or move randomly (if `randomize_tregs = 1`)
- Both use probabilistic movement based on local signal gradients

### 11. Add New Microbes (Leakage)
- **Pathogens:** Leak through injured epithelium (unless `sterile = 1`)
  - Rate: `rate_leak_pathogen_injury * mean(level_injury) * length(injury_site)`
  - Enter at y=1, weighted by injury level
- **Commensals:** Leak through epithelium (baseline + injury-enhanced)
  - Baseline: `rate_leak_commensal_baseline * grid_size`
  - Injury-enhanced: `rate_leak_commensal_injury * mean(level_injury) * length(injury_site)`

### 12. Update Phagocyte Phenotypes
- **M0 phagocytes:** Assess environment, potentially activate to M1 or M2
  - Calculate `danger_signal = avg_DAMPs + avg_PAMPs`
  - Calculate `avg_SAMPs`
  - Activate based on signal dominance
- **M1/M2 phagocytes:** After `active_age_limit` steps, reassess
  - With macrophage specificity: Use engulfment history
  - Without: Use environmental signals only

### 13. Update Treg Active Age
- Increment `active_age` for activated Tregs
- Deactivate Tregs after `active_age_limit` steps

### 14. Engulfment Process
- For each phagocyte:
  - Find colocated microbes
  - Attempt engulfment (probability = `activity_engulf`)
  - Update bacteria registry
  - Remove engulfed microbes from environment
  - Track deaths by phenotype

### 15. Treg Activation
- If `allow_tregs = 1`:
  - For each M1/M2 phagocyte with bacteria:
    - Find nearby Tregs (within `treg_vicinity_effect`)
    - Calculate commensal/pathogen ratio with discrimination
    - Activate Tregs if commensal ratio > threshold

### 16. Kill Microbes with ROS
- For each microbe:
  - Calculate average ROS in vicinity (`act_radius_ROS`)
  - Kill if ROS > `th_ROS_microbe`
  - Uses C++ function `kill_microbes_with_ros_cpp()` for 50-100x speedup

### 17. Update Epithelial Injury
- **Increase injury:**
  - From pathogens at epithelium: `+logistic_scaled(pathogen_count)`
  - From high ROS: `+1` if local ROS > `th_ROS_epith_injury`
- **Stochastic recovery:**
  - Each injured cell: probability `epith_recovery_chance` to reduce injury by 1
- **Cap:** Maximum injury level = 5

### 18. Save Abundances
- Record all population counts and cumulative death counters
- Append to longitudinal dataframe

---

## 8. Simulation Scenarios

The simulation explores multiple experimental conditions through scenario combinations:

### Scenario Parameters

1. **`control`:** (0 or 1)
   - 0: Normal simulation with ROS
   - 1: No ROS production (tests if infection resolves without ROS)

2. **`sterile`:** (0 or 1)
   - 0: Pathogens leak into lamina propria (infectious injury)
   - 1: No pathogen leakage (sterile injury)

3. **`allow_tregs`:** (0 or 1)
   - 0: Tregs present but inactive (cannot produce SAMPs)
   - 1: Tregs functional (can activate and produce SAMPs)

4. **`randomize_tregs`:** (0 or 1)
   - 0: Tregs chemotax toward DAMPs
   - 1: Tregs move randomly (no chemotaxis)

5. **`macspec_on`:** (0, 1, or 2)
   - 0: Vanilla model (no macrophage discrimination)
   - 1: Macrophages discriminate with same efficiency as Tregs
   - 2: Macrophages have perfect discrimination

### Scenario Filtering Rules

Not all combinations are scientifically meaningful. Invalid scenarios are filtered:

1. `allow_tregs = 0 AND randomize_tregs = 1` → INVALID
   - No point randomizing inactive Tregs

2. `macspec_on > 0 AND allow_tregs = 1 AND randomize_tregs = 1` → INVALID
   - Comparing macrophage discrimination requires proper Treg localization

3. `macspec_on > 0 AND allow_tregs = 1 AND randomize_tregs = 0` → INVALID
   - When comparing macrophage vs Treg discrimination, Tregs should be off

### Control Scenarios
Separate control scenarios with fixed parameters:
- `control = 1`
- `sterile = 0 or 1`
- `allow_tregs = 0`
- `randomize_tregs = 0`
- `macspec_on = 0`

### Total Simulations
For each parameter set:
- Valid scenarios: ~10-15 (after filtering)
- Replicates per scenario: 10
- Total simulations: ~100-150 per parameter set

---

## 9. Parameters

### Fixed Parameters (not in CSV)

Defined in `DLL_datagen_abm.R`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `t_max` | 5000 | Maximum simulation time steps |
| `grid_size` | 25 | Grid dimensions (25×25) |
| `num_reps` | 10 | Replicates per scenario |
| `n_phagocytes` | 125 | Number of phagocytes (20% of grid) |
| `n_tregs` | 125 | Number of Tregs (20% of grid) |
| `n_commensals_lp` | 20 | Initial commensals |
| `injury_percentage` | 60 | % of epithelium initially injured |
| `max_level_injury` | 5 | Maximum epithelial injury level |
| `max_cell_value_ROS` | 1 | Maximum ROS concentration |
| `max_cell_value_DAMPs` | 1 | Maximum DAMPs concentration |
| `max_cell_value_SAMPs` | 1 | Maximum SAMPs concentration |
| `max_cell_value_PAMPs` | 1 | Maximum PAMPs concentration |
| `act_radius_ROS` | 1 | Activation radius for ROS |
| `act_radius_treg` | 1 | Treg vicinity effect radius |
| `act_radius_DAMPs` | 1 | Activation radius for DAMPs |
| `act_radius_SAMPs` | 1 | Activation radius for SAMPs |
| `k_in` | 0.044 | Logistic function steepness |
| `x0_in` | 50 | Logistic function midpoint |
| `treg_vicinity_effect` | 1 | Distance for Treg-phagocyte interaction |
| `activity_ROS_M0_baseline` | 0 | M0 ROS production (none) |
| `activity_ROS_M2_baseline` | 0 | M2 ROS production (none) |

### Sampled Parameters (from CSV)

Read from `lhs_parameters_della.csv` (Latin Hypercube Sampling):

**Microbe leakage rates:**
- `rate_leak_pathogen_injury` - Pathogen leakage through injured epithelium
- `rate_leak_commensal_injury` - Commensal leakage through injured epithelium
- `rate_leak_commensal_baseline` - Baseline commensal leakage

**Epithelial parameters:**
- `epith_recovery_chance` - Probability of injury level reduction per time step

**Thresholds:**
- `th_ROS_microbe` - ROS concentration that kills microbes
- `th_ROS_epith_injury` - ROS concentration that injures epithelium
- `rat_com_pat_threshold` - Commensal ratio threshold for Treg activation
- `activation_threshold_danger` - Danger signal (DAMPs + PAMPs) threshold for M1 activation
- `activation_threshold_SAMPs` - SAMP threshold for M2 activation

**Diffusion speeds:**
- `diffusion_speed_DAMPs` - DAMPs diffusion rate
- `diffusion_speed_PAMPs` - PAMPs diffusion rate (NEW)
- `diffusion_speed_SAMPs` - SAMPs diffusion rate
- `diffusion_speed_ROS` - ROS diffusion rate

**Signal production:**
- `add_ROS` - ROS production rate by M1 phagocytes
- `add_DAMPs` - DAMPs production rate by injured epithelium
- `add_SAMPs` - SAMPs production rate by activated Tregs
- `add_PAMPs` - PAMPs production rate by pathogens (NEW)

**Decay rates:**
- `ros_decay` - ROS decay rate per time step
- `DAMPs_decay` - DAMPs decay rate per time step
- `SAMPs_decay` - SAMPs decay rate per time step
- `PAMPs_decay` - PAMPs decay rate per time step (NEW)

**Phagocyte activity:**
- `activity_engulf_M0_baseline` - M0 baseline engulfment probability (typically 0.05)
- `activity_engulf_M1_baseline` - M1 baseline engulfment probability
- `activity_engulf_M2_baseline` - M2 baseline engulfment probability
- `activity_ROS_M1_baseline` - M1 ROS production activity

**Phagocyte dynamics:**
- `cc_phagocyte` - Phagocyte carrying capacity (bacteria registry size)
- `active_age_limit` - Time steps before reassessing environment

**Discrimination:**
- `treg_discrimination_efficiency` - Treg ability to distinguish commensals from pathogens (0-1)
- `mac_discrimination_efficiency` - Derived from scenario and Treg discrimination (when `macspec_on > 0`)

**Recruitment:**
- `recruitment_rate_danger` - Macrophage recruitment rate from borders proportional to danger signal (0-0.5) (NEW)

---

## 10. Metrics Recorded (Longitudinal Data)

At each time step, the following data is recorded and saved in longitudinal dataframes:

### Epithelium
- **Injury distribution:** Number of cells at each injury level (0-5)
  - `epithelial_healthy` (level 0)
  - `epithelial_inj_1` through `epithelial_inj_5`

- **Epithelial health score:** Weighted sum of injury states
  - Formula: `epithelial_score = 6*healthy + 5*inj_1 + 4*inj_2 + 3*inj_3 + 2*inj_4 + 1*inj_5`
  - Higher score = healthier epithelium
  - Range: 25 (all maximally injured) to 150 (all healthy)

### Phagocytes
- `phagocyte_M0` - Number of resting M0 phagocytes
- `phagocyte_M1` - Number of pro-inflammatory M1 phagocytes
- `phagocyte_M2` - Number of anti-inflammatory M2 phagocytes

### Microbes
- `commensal` - Total number of live commensals
- `pathogen` - Total number of live pathogens

### Tregs
- `treg_resting` - Number of resting Tregs (phenotype 0)
- `treg_active` - Number of activated Tregs (phenotype 1)

### Microbe Death Counters (Cumulative)
Cumulative counts of microbes killed throughout simulation:

**Commensals killed:**
- `C_ROS` - By ROS
- `C_M0` - By M0 phagocytes
- `C_M1` - By M1 phagocytes
- `C_M2` - By M2 phagocytes

**Pathogens killed:**
- `P_ROS` - By ROS
- `P_M0` - By M0 phagocytes
- `P_M1` - By M1 phagocytes
- `P_M2` - By M2 phagocytes

### Steady State Detection
- **`time_ss`:** Time step when epithelial score reaches steady state
  - Calculated using `steady_state_idx()` function from `DATA_READ_FUNCTIONS.R`
  - Based on rolling window analysis (mean, SD, slope) of epithelial score
  - Uses 25% tail of trajectory to estimate asymptotic value

### Metadata
Each longitudinal dataframe includes:
- `t` - Time step
- `control` - Control scenario flag
- `sterile` - Sterile injury flag
- `tregs_on` - Tregs functional flag
- `macspec_on` - Macrophage specificity mode (0/1/2)
- `randomize_tregs` - Random Treg movement flag
- `param_set_id` - Parameter set identifier
- `rep_id` - Replicate number (0-9)

---

## 11. Implementation & Computational Optimization

### Script Organization

**Main script:** `DLL_datagen_abm.R`
- Reads parameters from CSV
- Manages scenario loops
- Calls simulation scripts

**Key MISC scripts:**
- `ASSIGN_PARAMETERS.R` - Assigns parameters from CSV to simulation variables
- `RUN_REPS_CPP_ABM_PAMPS.R` - Main simulation loop with PAMPs support (CURRENT VERSION)
- `RUN_REPS_CPP_ABM.R` - Legacy simulation loop (without PAMPs)
- `FAST_FUNCTIONS_CPP.R` - Loads C++ accelerated functions
- `FAST_FUNCTIONS.cpp` - C++ implementations of computationally intensive operations
- `PLOT_FUNCTIONS_ABM.R` - Visualization functions
- `DATA_READ_FUNCTIONS.R` - Data processing and steady-state detection
- `CONVERT_TO_DATAFRAME_ABM.R` - Converts simulation state to dataframes for plotting

### C++ Acceleration

**Motivation:** Pure R implementation is too slow for parameter exploration

**Implementation:** Using Rcpp to compile C++ code that interfaces with R

**Accelerated functions:** (20-100x speedup)

| Function | Speedup | Purpose |
|----------|---------|---------|
| `diffuse_matrix_cpp()` | 5-10x | Matrix diffusion (Laplacian) |
| `kill_microbes_with_ros_cpp()` | 50-100x | Batch microbe death calculation |
| `calculate_phagocyte_signals_cpp()` | 20-50x | Batch signal averaging for phagocytes |
| `calculate_epithelial_ros_cpp()` | 10-20x | Batch ROS calculation for epithelium |
| `update_SAMPs_batch_cpp()` | 20-50x | Batch SAMP production by Tregs |
| `shift_insert_fast_cpp()` | 5-10x | Bacteria registry updates |
| `find_nearby_tregs_cpp()` | 10-20x | Spatial proximity queries |

**Fallback:** R implementations available if C++ compilation fails

**Status checking:** Use `check_cpp_status()` to see which functions are accelerated

### Parallelization

**Job parallelization:** `DLL_datagen_abm.R` accepts command-line arguments
- Argument 1: Total number of chunks
- Argument 2: Chunk index to process
- Allows parallel execution across computing cluster

**Example:**
```bash
Rscript DLL_datagen_abm.R 10 1  # Process chunk 1 of 10
Rscript DLL_datagen_abm.R 10 2  # Process chunk 2 of 10
...
```

### Output
- Results saved as RDS files: `longitudinal_df_param_set_id_{ID}_control_{0/1}_sterile_{0/1}_macspec_{0/1/2}_tregs_{0/1}_trnd_{0/1}.rds`
- Output directory: `/scratch/gpfs/CMETCALF/sim_abm` (configurable in script)

---

## 12. Key Changes from Original Documentation

### Major Additions

1. **PAMPs (Pathogen-Associated Molecular Patterns):**
   - Completely new signaling molecule
   - Released by pathogens
   - Combines with DAMPs to form danger signal
   - Has own diffusion, decay, and production parameters
   - Enables distinction between pathogenic infection and sterile injury

2. **Macrophage Specificity (`macspec_on`):**
   - New feature allowing macrophages to discriminate between commensals and pathogens
   - Three modes: off (0), equal to Tregs (1), perfect (2)
   - Uses engulfment history for polarization decisions
   - Implements danger-biased logic (M1 easier to trigger than M2)

3. **C++ Acceleration:**
   - Entire computational pipeline optimized with C++ implementations
   - 20-100x speedup for most intensive operations
   - Essential for large parameter explorations

4. **Macrophage Recruitment from Borders:**
   - Completely new feature for dynamic immune cell recruitment
   - Location-specific recruitment proportional to local danger signals (DAMPs + PAMPs)
   - Macrophages enter from top, left, and right borders (not epithelium)
   - Poisson-distributed stochastic recruitment
   - Allows model to capture dynamic changes in macrophage numbers during inflammation
   - Controlled by `recruitment_rate_danger` parameter (range: 0-0.5)

### Logic Changes

1. **Danger Signal Definition:**
   - Old: Just DAMPs
   - New: `danger_signal = DAMPs + PAMPs`
   - Affects phagocyte activation and chemotaxis

2. **Phagocyte Activation Threshold:**
   - Old: `activation_threshold_DAMPs`
   - New: `activation_threshold_danger` (for combined DAMPs + PAMPs)

3. **Microbe Movement Order:**
   - Now occurs BEFORE DAMP calculation
   - Ensures DAMPs reflect current microbe positions

4. **DAMP Generation from Microbes:**
   - Now explicitly includes both commensals and pathogens at epithelium
   - Uses `logistic_scaled_0_to_5_quantized()` for scaling

5. **Phagocyte Chemotaxis:**
   - Old: Toward DAMPs only
   - New: Toward `DAMPs + PAMPs`

### Scenario Complexity

1. **Control Scenarios:** Explicitly separated (no ROS to test necessity)
2. **Macrophage Specificity:** New dimension for comparison with Treg specificity
3. **Scenario Filtering:** More sophisticated rules to avoid meaningless combinations

### Data Output

1. **Epithelial Health Score:** New aggregate metric for epithelial state
2. **Steady State Detection:** Automated detection of when system reaches equilibrium
3. **Enhanced Metadata:** Tracks all scenario parameters in longitudinal data

---

## 13. Biological Interpretation & Model Logic

### Immune Response Phases

The model captures the progression of immune responses:

**Phase 1: Injury & Danger Signals (t = 0-100)**
- Epithelium injured
- DAMPs released from injured cells
- Microbes (commensals + pathogens) leak through damaged barrier
- PAMPs released by pathogens (if present)

**Phase 2: Phagocyte Activation (t = 50-500)**
- Phagocytes chemotax toward danger (DAMPs + PAMPs)
- M0 → M1 activation in high danger zones
- M1 phagocytes produce ROS
- ROS kills microbes AND causes collateral epithelial damage

**Phase 3: Treg Response (t = 100-1000)**
- Tregs chemotax toward DAMPs
- Encounter M1/M2 phagocytes presenting antigens
- If commensals dominate: Tregs activate
- Activated Tregs produce SAMPs

**Phase 4: Resolution (t = 500-5000)**
- SAMPs shift balance toward M2
- M2 phagocytes clear microbes without ROS
- Epithelium recovers stochastically
- System reaches homeostasis

### Model Predictions

**Treg Benefits (when helping):**
- Reduce collateral damage by suppressing excessive M1/ROS
- Promote M2-mediated clearance
- Allow commensal tolerance
- Accelerate epithelial recovery

**Treg Costs (when harming):**
- Suppress necessary M1 responses to pathogens
- Allow pathogen persistence
- Delay clearance of infection

**Macrophage Specificity Benefits:**
- Can distinguish commensals from pathogens based on engulfment
- M2 polarization only with concordant safety signals
- May reduce need for Tregs in some contexts

### Parameter Space Exploration

The Latin Hypercube Sampling (LHS) of parameters explores:
- Different kinetics of injury, recovery, and clearance
- Trade-offs between antimicrobial efficacy and collateral damage
- Scenarios where Tregs help, harm, or are irrelevant
- Role of macrophage discrimination

---

## Appendix: Function Reference

### Helper Functions (from FAST_FUNCTIONS_CPP.R)

- **`get_middle_percent(seq_vector, percent)`:** Returns middle X% of a sequence
  - Used to define injury site

- **`logistic_scaled_0_to_5_quantized(x, k, x0)`:** Maps counts to injury levels (0-5)
  - Uses logistic function for smooth scaling
  - Quantized to integer levels

- **`iszero_coordinates(x)`:** Generates orthogonal movement when primary direction is 0
  - Ensures microbes don't stay stationary at epithelium

- **`sample_rbeta(alpha, beta)`:** Beta distribution sampling
  - Used in discrimination calculations (commented out in current version)

- **`steady_state_idx(x, k, tail_frac, tol_abs, tol_sd, tol_slope)`:** Detects steady state
  - Rolling window analysis
  - Checks: mean near asymptote, low SD, near-zero slope
  - From `DATA_READ_FUNCTIONS.R`

### C++ Functions (from FAST_FUNCTIONS.cpp)

All functions are vectorized for batch operations where possible. Fallback R implementations available in FAST_FUNCTIONS_CPP.R if compilation fails.

---

## Document Metadata

**Version:** 2.1 (Updated)
**Date:** 2025-12-16
**Primary Changes:** Added PAMPs system, macrophage specificity, C++ acceleration, border recruitment system, complete simulation algorithm documentation, updated parameter names, clarified logic
**Corresponds to:** `DLL_datagen_abm.R` and `MISC/RUN_REPS_CPP_ABM_PAMPS.R`
