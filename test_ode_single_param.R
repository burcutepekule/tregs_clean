#!/usr/bin/env Rscript
# Test script to run ODE system with a single parameter set
# This validates the ODE implementation before running full data generation

rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(deSolve)
library(ggplot2)

cat("=========================================\n")
cat("ODE SYSTEM TEST - SINGLE PARAMETER SET\n")
cat("=========================================\n\n")

source("./MISC/ODE_SYSTEM.R")
source("./MISC/DATA_READ_FUNCTIONS.R")
source("./MISC/FAST_FUNCTIONS_CPP.R")

# ============================================================================
# LOAD PARAMETERS
# ============================================================================

cat("Loading parameters...\n")
params_df = read.csv("./lhs_parameters_della.csv", stringsAsFactors = FALSE)
cat("Total parameter sets available:", nrow(params_df), "\n\n")

# Use first parameter set for testing
param_set_use = params_df %>% dplyr::filter(param_set_id == 0)

# ============================================================================
# FIXED PARAMETERS
# ============================================================================

t_max      = 5000
num_reps   = 3  # Just 3 replicates for testing

grid_size       = 25
n_phagocytes    = round(grid_size*grid_size*0.20)
n_tregs         = round(grid_size*grid_size*0.20)
n_commensals_lp = 20

injury_percentage = 60
max_level_injury  = 5

k_in  = 0.044
x0_in = 50

# ============================================================================
# SCENARIO SETUP
# ============================================================================

sterile         = 0
allow_tregs     = 0
randomize_tregs = 0
macspec_on      = 0
ros_level       = 0
pat_level       = 2  # Medium pathogen level

cat("Test scenario:\n")
cat("  Parameter set ID:", param_set_use$param_set_id, "\n")
cat("  Pathogen level:", pat_level, "\n")
cat("  Tregs:", allow_tregs, "\n")
cat("  Macrophage specificity:", macspec_on, "\n\n")

# ============================================================================
# ASSIGN PARAMETERS
# ============================================================================

cat("Assigning parameters...\n")
source("./MISC/ASSIGN_PARAMETERS.R")

# ============================================================================
# PREPARE ODE PARAMETERS
# ============================================================================

ode_params <- list(
  # === Grid/environment parameters ===
  grid_size = grid_size,
  E_max = grid_size^2 * 6,
  n_phagocytes = n_phagocytes,
  n_commensals_lp = n_commensals_lp,
  injury_percentage = injury_percentage,
  max_level_injury = max_level_injury,

  # === Macrophage polarization ===
  activation_threshold_danger = activation_threshold_danger,
  activation_threshold_SAMPs = activation_threshold_SAMPs,
  polarization_speed = 0.1,
  K_polarization = 0.5,

  # === Macrophage recruitment ===
  recruitment_rate_danger = recruitment_rate_danger,

  # === Macrophage activities ===
  activity_ROS_M0_baseline = activity_ROS_M0_baseline,
  activity_ROS_M1_baseline = activity_ROS_M1_baseline,
  activity_ROS_M2_baseline = activity_ROS_M2_baseline,
  activity_engulf_M0_baseline = activity_engulf_M0_baseline,
  activity_engulf_M1_baseline = activity_engulf_M1_baseline,
  activity_engulf_M2_baseline = activity_engulf_M2_baseline,

  # === Epithelial dynamics ===
  epith_recovery_chance = epith_recovery_chance,
  damage_rate_pathogen = 1.0,
  damage_rate_ROS = 2.0,
  K_pathogen_damage = 50,
  th_ROS_epith_injury = th_ROS_epith_injury,

  # === Microbe leakage ===
  rate_leak_pathogen_injury = rate_leak_pathogen_injury,
  rate_leak_commensal_baseline = rate_leak_commensal_baseline,
  rate_leak_commensal_injury = rate_leak_commensal_injury,

  # === Microbe killing ===
  th_ROS_microbe = th_ROS_microbe,
  kill_rate_ROS = 2.0,
  K_phagocytosis = 100,

  # === Treg parameters ===
  allow_tregs = allow_tregs,
  treg_activation_efficiency = 0.01,
  K_treg_activation = 50,

  # === Signaling molecule production ===
  add_ROS = add_ROS,
  add_DAMPs = add_DAMPs,
  add_PAMPs = add_PAMPs,
  add_SAMPs = add_SAMPs,

  # === Signaling molecule decay ===
  ros_decay = ros_decay,
  DAMPs_decay = DAMPs_decay,
  PAMPs_decay = PAMPs_decay,
  SAMPs_decay = SAMPs_decay,

  # === Half-saturation constants ===
  K_DAMPS_microbe = 100
)

cat("ODE parameters prepared\n\n")

# ============================================================================
# RUN ODE SIMULATIONS
# ============================================================================

cat("Running ODE simulations...\n")
longitudinal_df_keep = data.frame()

for (reps_in in 0:(num_reps-1)) {
  cat(sprintf("  Replicate %d/%d...", reps_in + 1, num_reps))

  # Run ODE
  start_time = Sys.time()
  ode_output <- run_ode_simulation(
    parameters = ode_params,
    t_max = t_max,
    dt = 1
  )

  # Process output
  scenario_info <- list(
    param_set_id = param_set_use$param_set_id,
    rep_id = reps_in,
    sterile = sterile,
    macspec_on = macspec_on,
    allow_tregs = allow_tregs,
    ros_level = ros_level,
    pat_level = pat_level,
    randomize_tregs = randomize_tregs
  )

  longitudinal_df <- process_ode_output(ode_output, ode_params, scenario_info)
  longitudinal_df_keep <- rbind(longitudinal_df_keep, longitudinal_df)

  end_time = Sys.time()
  elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(" %.2f seconds ✓\n", elapsed))
}

cat("\nSimulations complete!\n\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("=========================================\n")
cat("SUMMARY STATISTICS\n")
cat("=========================================\n\n")

# Get final timepoint (average across replicates)
final_df <- longitudinal_df_keep %>%
  filter(t == t_max) %>%
  summarise(
    epithelial_score = mean(epithelial_score),
    pathogen = mean(pathogen),
    commensal = mean(commensal),
    phagocyte_M0 = mean(phagocyte_M0),
    phagocyte_M1 = mean(phagocyte_M1),
    phagocyte_M2 = mean(phagocyte_M2),
    ROS = mean(ROS),
    DAMPS = mean(DAMPS),
    PAMPS = mean(PAMPS),
    SAMPS = mean(SAMPS)
  )

cat("Final state (t =", t_max, ", averaged across replicates):\n")
cat(sprintf("  Epithelial score: %.2f / %.0f (%.1f%% healthy)\n",
            final_df$epithelial_score,
            ode_params$E_max,
            100 * final_df$epithelial_score / ode_params$E_max))
cat(sprintf("  Pathogens: %.2f\n", final_df$pathogen))
cat(sprintf("  Commensals: %.2f\n", final_df$commensal))
cat(sprintf("  Macrophages: M0=%.1f, M1=%.1f, M2=%.1f\n",
            final_df$phagocyte_M0, final_df$phagocyte_M1, final_df$phagocyte_M2))
cat(sprintf("  Signals: ROS=%.3f, DAMPS=%.3f, PAMPS=%.3f, SAMPS=%.3f\n",
            final_df$ROS, final_df$DAMPS, final_df$PAMPS, final_df$SAMPS))

# Time to steady state
ss_df <- longitudinal_df_keep %>%
  group_by(rep_id) %>%
  summarise(
    time_ss_e = first(time_ss_e),
    time_ss_p = first(time_ss_p)
  )

cat("\nSteady state times:\n")
cat(sprintf("  Epithelial: %.0f ± %.0f timesteps\n",
            mean(ss_df$time_ss_e, na.rm = TRUE),
            sd(ss_df$time_ss_e, na.rm = TRUE)))
cat(sprintf("  Pathogen: %.0f ± %.0f timesteps\n",
            mean(ss_df$time_ss_p, na.rm = TRUE),
            sd(ss_df$time_ss_p, na.rm = TRUE)))

# Peak pathogen load
peak_df <- longitudinal_df_keep %>%
  group_by(rep_id) %>%
  summarise(peak_pathogen = max(pathogen))

cat("\nPeak pathogen load:\n")
cat(sprintf("  %.2f ± %.2f pathogens\n",
            mean(peak_df$peak_pathogen),
            sd(peak_df$peak_pathogen)))

cat("\n=========================================\n")
cat("TEST COMPLETED SUCCESSFULLY!\n")
cat("=========================================\n\n")

cat("Output saved to: test_ode_output.rds\n")
saveRDS(longitudinal_df_keep, "test_ode_output.rds")

# ============================================================================
# OPTIONAL: CREATE SIMPLE PLOT
# ============================================================================

cat("\nCreating diagnostic plots...\n")

# Average across replicates for cleaner visualization
avg_df <- longitudinal_df_keep %>%
  group_by(t) %>%
  summarise(
    epithelial_score = mean(epithelial_score),
    pathogen = mean(pathogen),
    commensal = mean(commensal),
    M0 = mean(phagocyte_M0),
    M1 = mean(phagocyte_M1),
    M2 = mean(phagocyte_M2),
    ROS = mean(ROS),
    DAMPS = mean(DAMPS),
    PAMPS = mean(PAMPS),
    SAMPS = mean(SAMPS)
  )

# Plot 1: Populations
p1 <- ggplot(avg_df, aes(x = t)) +
  geom_line(aes(y = pathogen, color = "Pathogens"), linewidth = 1) +
  geom_line(aes(y = commensal, color = "Commensals"), linewidth = 1) +
  geom_line(aes(y = M1, color = "M1 macrophages"), linewidth = 1) +
  geom_line(aes(y = M2, color = "M2 macrophages"), linewidth = 1) +
  labs(title = "Population Dynamics",
       x = "Time (steps)",
       y = "Count",
       color = "Population") +
  theme_minimal()

# Plot 2: Epithelial health
p2 <- ggplot(avg_df, aes(x = t, y = epithelial_score)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_hline(yintercept = ode_params$E_max, linetype = "dashed", color = "gray50") +
  labs(title = "Epithelial Health Score",
       x = "Time (steps)",
       y = "Health Score") +
  theme_minimal()

# Plot 3: Signaling molecules
p3 <- ggplot(avg_df, aes(x = t)) +
  geom_line(aes(y = ROS, color = "ROS"), linewidth = 1) +
  geom_line(aes(y = DAMPS, color = "DAMPs"), linewidth = 1) +
  geom_line(aes(y = PAMPS, color = "PAMPs"), linewidth = 1) +
  geom_line(aes(y = SAMPS, color = "SAMPs"), linewidth = 1) +
  labs(title = "Signaling Molecules",
       x = "Time (steps)",
       y = "Concentration",
       color = "Signal") +
  theme_minimal()

ggsave("test_ode_populations.png", p1, width = 10, height = 6, dpi = 300)
ggsave("test_ode_epithelial.png", p2, width = 10, height = 6, dpi = 300)
ggsave("test_ode_signals.png", p3, width = 10, height = 6, dpi = 300)

cat("Plots saved:\n")
cat("  - test_ode_populations.png\n")
cat("  - test_ode_epithelial.png\n")
cat("  - test_ode_signals.png\n\n")

cat("Test complete!\n")
