# CMEE 2024 HPC exercises R code pro forma
# For neutral model cluster run

rm(list=ls()) # good practice 
source("abc123_HPC_2024_main.R")

# CMEE 2024 HPC Cluster Run Script for Neutral Model Simulations

# Read job number from cluster
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#iter <- 1

# Set a unique random seed for each job
set.seed(iter)

# Assign community sizes based on job number
if (iter <= 25) {
  size <- 500
} else if (iter <= 50) {
  size <- 1000
} else if (iter <= 75) {
  size <- 2500
} else {
  size <- 5000
}

# Define other simulation parameters
speciation_rate <- 0.1  
wall_time <- 690   
interval_rich <- 1       
interval_oct <- size / 10  
burn_in_generations <- 8 * size  

# Generate unique output filename
output_file_name <- paste0("neutral_cluster_output_", iter, ".rda")

# Run the neutral simulation
neutral_cluster_run(speciation_rate, size, wall_time, interval_rich, 
                    interval_oct, burn_in_generations, output_file_name)
