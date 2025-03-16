# Clear the workspace and turn off graphics
rm(list = ls())
graphics.off()

# Source the file that provides necessary functions (e.g., stochastic_simulation, state_initialise_adult/spread)
source("demographic.R")
source("yh4724_HPC_2024_main.R")

# Read the job number from the environment variable PBS_ARRAY_INDEX
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
if (is.na(iter)) {
  iter <- 1  # For local testing, if PBS_ARRAY_INDEX is not set
}

# Set simulation parameters
simulation_length <- 120  # Simulation length in time steps (e.g., 120 steps = 10 years if time scale is in months)
num_stages <- 4           # There are 4 life stages in the model

# Define the growth matrix and reproduction matrix (same as used in Question 2)
growth_matrix <- matrix(c(0.1, 0.0, 0.0, 0.0,
                          0.5, 0.4, 0.0, 0.0,
                          0.0, 0.4, 0.7, 0.0,
                          0.0, 0.0, 0.25, 0.4), nrow = 4, byrow = TRUE)

reproduction_matrix <- matrix(c(0.0, 0.0, 0.0, 2.6,
                                0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0), nrow = 4, byrow = TRUE)

# Clutch distribution as provided
clutch_distribution <- c(0.06, 0.08, 0.13, 0.15, 0.16, 0.18, 0.15, 0.06, 0.03)

# Determine the initial condition based on the job number iter
if (iter >= 1 && iter <= 25) {
  init <- state_initialise_adult(num_stages, 100)
  init_label <- "big_adult"
} else if (iter >= 26 && iter <= 50) {
  init <- state_initialise_adult(num_stages, 10)
  init_label <- "small_adult"
} else if (iter >= 51 && iter <= 75) {
  init <- state_initialise_spread(num_stages, 100)
  init_label <- "big_spread"
} else {
  init <- state_initialise_spread(num_stages, 10)
  init_label <- "small_spread"
}

# Set the random number seed so that each parallel simulation is unique
set.seed(iter)

# Create a filename to store simulation results; appending iter ensures unique filenames
outfile <- paste0("results_", iter, ".rda")

# Initialize a list to store the results of 150 simulations
results_list <- list()

# Run 150 simulations using stochastic_simulation and store each result in results_list
for(i in 1:150) {
  sim_result <- stochastic_simulation(init, growth_matrix, reproduction_matrix, clutch_distribution, simulation_length)
  results_list[[i]] <- sim_result
}

# Save the results list to the specified output file
save(results_list, file = outfile)

