# ============================================================================
# MAIN SIMULATION LOOP-DIFFUSION-BASED MODEL
# ============================================================================
# This replaces the agent-based model with a density-field approach where
# microbes and lymphocytes are represented as continuous density grids
# that evolve via diffusion, decay, and reaction terms.
#
# Key simplifications:
#-No individual agent tracking (no bacteria registry, no agent IDs)
#-All entities represented as density fields on the same 25x25 grid
#-Phenotype updates happen every active_age_limit steps
#-Much faster execution while preserving spatial dynamics
# ============================================================================

# ============================================================================
# MAIN SIMULATION LOOP - DIFFUSION MODEL
# ============================================================================
# for (t in 1:delay_response) { # DEPRICIATED
#   source('./MISC/MAIN_SIM_DIFFUSION_PRE_IMMUNE.R')
# }

for (t in (delay_response+1):t_max) {
  source('./MISC/MAIN_SIM_DIFFUSION_RPAT.R')
}
