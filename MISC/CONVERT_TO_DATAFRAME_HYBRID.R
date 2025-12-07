# Create the complete dataframe with bacteria registry
phagocytes <- data.frame(
  x = phagocyte_x,
  y = phagocyte_y,
  pathogens_engulfed = phagocyte_pathogen_memory,#engulfed to memory
  commensals_engulfed = phagocyte_commensal_memory,
  phenotype = phagocyte_phenotype,
  activity_ROS = phagocyte_activity_ROS,
  activity_engulf = phagocyte_activity_engulf,
  active_age = phagocyte_active_age
)

# Add the bacteria registry as a list column
phagocytes$bacteria_registry <- 1

tregs=data.frame(
  x = treg_x,
  y = treg_y,
  active_age = treg_active_age,
  phenotype = treg_phenotype,
  activity_SAMPs_binary = treg_activity_SAMPs_binary
)

