# ABM vs HYBRID Simulation Discrepancy Analysis

## Summary of Observed Discrepancy

**ABM Model (left panel, tregs_on=0):** Shows significant epithelial cell death and fluctuation (count drops from ~150 to ~90-110 with high variance)

**HYBRID Model (left panel, tregs_on=0):** Shows stable epithelial counts (~150) with minimal fluctuation

## Root Causes Identified

### 1. **CRITICAL BUG: Coordinate System Mismatch in HYBRID Model**

#### ABM Coordinate System (Correct):
- Epithelium: y=0 (conceptual), first tissue layer at y=1
- Microbes added at: **y=1** (at epithelium)
- Immune cells can access: y≥1
- Epithelial damage calculated from pathogens at: **y=1**
- DAMPs from microbes calculated at: **y=1**

**Relevant code (RUN_REPS_CPP_ABM.R):**
```r
Line 19:  epithelium y = rep(0, grid_size)  # Conceptual position
Line 50:  pathogen_coords at y=1            # At epithelium
Line 137: epithelium_pathogens = pathogen_coords[, "y"] == 1  # Counts at y=1
Line 155: DAMPs[1, ] += epithelium$level_injury * add_DAMPs  # y=1
Line 158: DAMPs from microbes at y=1
Line 622: Injury from pathogen_epithelium_counts (from y=1)
```

#### HYBRID Coordinate System (INCORRECT):
- Epithelium: **y=1**
- Microbes added at: **y=2** ❌ (one layer ABOVE epithelium!)
- Immune cells can access: y≥2
- Epithelial damage calculated from pathogens at: **y=2** ❌ (wrong layer!)
- DAMPs from microbes calculated at: **y=2** ❌ (wrong layer!)

**Buggy code (RUN_REPS_CPP_HYBRID.R):**
```r
Line 29:  epithelium y = rep(1, grid_size)  # Epithelium at y=1
Line 60:  P_field[2, x_start] += 1          # ❌ Pathogens at y=2, not y=1!
Line 116: microbe_load_epithelium = P_field[2, ] + C_field[2, ]  # ❌ Wrong row!
Line 130: P_field[2, x] += n_new_pathogens  # ❌ Adding at y=2, not y=1!
Line 482: pathogen_load = P_field[2, x]     # ❌ Checking y=2, not y=1!
```

### 2. **Impact of the Bug**

The HYBRID model adds pathogens at y=2 (one layer above the epithelium at y=1), causing:

1. **Reduced epithelial damage**: Pathogens at y=2 don't directly contact epithelium at y=1
2. **Reduced DAMP production**: Microbe-induced DAMPs calculated from wrong layer
3. **Stable epithelial counts**: Without proper pathogen damage, epithelium remains healthy
4. **Reduced immune response**: Lower DAMPs → less macrophage activation → less pathogen clearing

This creates a vicious cycle where the system appears stable because pathogens aren't damaging the epithelium they should be touching.

### 3. **Microbe Diffusion Does Not Compensate**

Even though microbes diffuse with `D_microbe = 0.001` and `reflect_top=TRUE`:
- Diffusion coefficient is very small (0.001)
- Reflective boundary at y=1 creates no-flux condition
- Microbes would need many timesteps to diffuse from y=2 to y=1
- Even with diffusion, the **damage calculation still checks the wrong row** (line 482: `P_field[2, x]`)

### 4. **Additional Coordinate-Related Issues**

#### Movement Boundaries:
**ABM:**
```r
Line 208: y_range = max(1, y - 1):min(grid_size, y + 1)  # Tregs can go to y=1
Line 238: y_range = max(1, y - 1):min(grid_size, y + 1)  # Phagocytes can go to y=1
```

**HYBRID:**
```r
Line 163: y_range = max(2, y - 1):min(grid_size, y + 1)  # Tregs blocked from y=1
Line 197: y_range = max(2, y - 1):min(grid_size, y + 1)  # Phagocytes blocked from y=1
```

This prevents immune cells from reaching the epithelium layer in HYBRID.

### 5. **Memory Decay Implementation Difference**

**ABM (discrete decay):**
```r
Line 495: phagocyte_pathogens_engulfed[i] = max(0, ...-1*(t%%3==0))  # Decrease by 1 every 3 steps
```

**HYBRID (exponential decay):**
```r
Line 228: phagocyte_pathogen_memory = phagocyte_pathogen_memory * memory_decay_coeff
# memory_decay_coeff = 0.01 means memory drops to 1% per timestep!
```

This is an extremely fast decay (99% loss per timestep), which may cause macrophages to "forget" what they engulfed almost immediately.

## Required Fixes for HYBRID Model

### Fix 1: Correct Microbe Addition Location
```r
# Line 60 (initial pathogens):
P_field[1, x_start] = P_field[1, x_start] + 1  # Changed from [2,] to [1,]

# Line 68 (initial commensals):
y_start = sample(1:grid_size, 1)  # Changed from 2:grid_size

# Lines 130, 134 (continuous entry):
P_field[1, x] = P_field[1, x] + n_new_pathogens  # Changed from [2,]
C_field[1, x] = C_field[1, x] + n_new_commensals  # Changed from [2,]
```

### Fix 2: Correct DAMP Calculation
```r
# Line 113 (DAMPs at epithelium):
DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury * add_DAMPs  # Changed from [2,]

# Line 116 (DAMPs from microbe contact):
microbe_load_epithelium = P_field[1, ] + C_field[1, ]  # Changed from [2,]
DAMPs[1, ] = DAMPs[1, ] + logistic_scaled_0_to_5_quantized(microbe_load_epithelium) * add_DAMPs
```

### Fix 3: Correct Epithelial Injury Calculation
```r
# Line 482 (injury from pathogens):
pathogen_load = P_field[1, x]  # Changed from P_field[2, x]
```

### Fix 4: Allow Immune Cells to Access Epithelium Layer
```r
# Lines 163, 197 (movement range):
y_range = max(1, y - 1):min(grid_size, y + 1)  # Changed from max(2,...)
```

### Fix 5: Adjust Memory Decay (Optional)
```r
# Line 98 (in datagen file):
memory_decay_coeff = 0.99  # Or use discrete decay like ABM: -1*(t%%3==0)
```

## Verification Steps

After implementing fixes:
1. Run both models with same parameters
2. Compare epithelial dynamics - should now show similar fluctuation patterns
3. Check pathogen/commensal counts - should be comparable
4. Verify macrophage activation patterns align

## Expected Outcome

With these fixes, the HYBRID model should:
- Show epithelial cell death when tregs_on=0 (similar to ABM)
- Produce comparable DAMP/SAMP dynamics
- Have pathogens actually damaging the epithelium
- Maintain computational efficiency while matching ABM results
