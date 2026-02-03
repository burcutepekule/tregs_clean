plot_on    = 0
plot_every = 0
t_max      = 2500
num_reps   = 10

grid_size       = 25
n_phagocytes    = round(grid_size*grid_size*0.20)
n_tregs         = round(grid_size*grid_size*0.20)
n_commensals_lp = 20
max_total_phagocytes = round(grid_size*grid_size*0.80)

injury_percentage = 60

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1
max_cell_value_PAMPs = 1

act_radius_ROS   = 1
act_radius_treg  = 1
act_radius_DAMPs = 1
act_radius_SAMPs = 1
act_radius_PAMPs = 1

# Logistic function parameters (for epithelial injury calculation)
### ORIGINAL PARAMETERS
k_in_log  = 0.044
x0_in_log = 50
c_in_log  = 5
# max_level_injury     = 5
# initial_injury_level = 1

# # Exponential function parameters (for epithelial injury calculation)
k_in  = 0.01
x0_in = 5
max_level_injury     = 2*x0_in # that's where both functions are converging from pathogen and ROS based injury! - this is per cell
initial_injury_level = 1

# Bacteria natural lifespan (in time steps)
age_max_bacteria = 100