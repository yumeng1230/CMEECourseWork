# CMEE 2024 HPC exercises R code main pro forma
# You don't HAVE to use this but it will be very helpful.
# If you opt to write everything yourself from scratch please ensure you use
# EXACTLY the same function and parameter names and beware that you may lose
# marks if it doesn't work properly because of not using the pro-forma.

name <- "Yumeng Huang"
preferred_name <- "Yumeng"
email <- "yh4724@ic.ac.uk"
username <- "yh4724"

# Please remember *not* to clear the work space here, or anywhere in this file.
# If you do, it'll wipe out your username information that you entered just
# above, and when you use this file as a 'toolbox' as intended it'll also wipe
# away everything you're doing outside of the toolbox.  For example, it would
# wipe away any automarking code that may be running and that would be annoying!

library(ggplot2)
# Section One: Stochastic demographic population model

# Question 0

state_initialise_adult <- function(num_stages, initial_size) {
  state_vector <- rep(0, num_stages)  # Create a vector of zeros with length num_stages
  state_vector[num_stages] <- initial_size  # Assign all individuals to the last stage (adult)
  return(state_vector)
}

state_initialise_spread <- function(num_stages, initial_size) {
  base_number <- floor(initial_size / num_stages)  # Base number of individuals per stage
  remainder <- initial_size %% num_stages  # Remaining individuals to be distributed
  state_vector <- rep(base_number, num_stages)  # Create a vector with base numbers
  
  if (remainder > 0) {
    state_vector[1:remainder] <- state_vector[1:remainder] + 1  # Distribute the remainder starting from the youngest stage
  }
  
  return(state_vector)
}

sum_vect <- function(a, b) {
  if (length(a) != length(b)) {
    stop("Vectors must be the same length")
  }
  return(a + b)
}

# Question 1
question_1 <- function(){
  
  # Load the demographic model functions
  source("Demographic.R")
  
  # Define the projection matrix
  growth_matrix <- matrix(c(0.1, 0.0, 0.0, 0.0,
                            0.5, 0.4, 0.0, 0.0,
                            0.0, 0.4, 0.7, 0.0,
                            0.0, 0.0, 0.25, 0.4),
                          nrow = 4, ncol = 4, byrow = TRUE)
  
  reproduction_matrix <- matrix(c(0.0, 0.0, 0.0, 2.6,
                                  0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0),
                                nrow = 4, ncol = 4, byrow = TRUE)
  
  # Compute the full projection matrix
  projection_matrix <- reproduction_matrix + growth_matrix
  
  # Function to simulate population projection over time
  project_population <- function(initial_state, projection_matrix, time_steps) {
    num_stages <- length(initial_state)
    population_over_time <- matrix(0, nrow = time_steps + 1, ncol = num_stages)
    population_over_time[1, ] <- initial_state  # Set initial population distribution
    
    for (t in 1:time_steps) {
      # Apply matrix multiplication to project population dynamics
      population_over_time[t + 1, ] <- projection_matrix %*% population_over_time[t, ]
    }
    
    return(rowSums(population_over_time))  # Return total population size over time
  }
  
  # Define simulation parameters
  num_stages <- 4
  time_steps <- 24
  initial_population_adult <- c(0, 0, 0, 100)  # 100 individuals in the adult stage
  initial_population_spread <- state_initialise_spread(num_stages, 100)  # 100 individuals spread across stages
  
  # Run population simulations
  pop_adult <- project_population(initial_population_adult, projection_matrix, time_steps)
  pop_spread <- project_population(initial_population_spread, projection_matrix, time_steps)
  
  # Create and save the population growth comparison plot
  png(filename="question_1.png", width = 600, height = 400)
  plot(0:time_steps, pop_adult, type = "l", col = "red", lwd = 2, ylim = range(c(pop_adult, pop_spread)),
       xlab = "Time (Years)", ylab = "Total Population Size", main = "Population Growth Comparison")
  lines(0:time_steps, pop_spread, col = "blue", lwd = 2)
  legend("topright", legend = c("All Adults", "Evenly Spread"), col = c("red", "blue"), lwd = 2)
  Sys.sleep(0.1)  
  dev.off()
  
  # Return explanation of the results
  return("The initial distribution of individuals across life stages affects both short-term and long-term population growth. 
  When all individuals start as adults, reproduction begins immediately, leading to a rapid initial increase in population. 
  In contrast, when individuals are spread across different life stages, there is a delay in reproduction as younger individuals need time to mature. 
  This results in slower initial growth, but eventually, the population stabilizes with a more balanced age structure.")
}


# Question 2
question_2 <- function(){
  
  # Load the required demographic model functions
  source("Demographic.R")
  
  # Define the growth and reproduction matrices separately
  growth_matrix <- matrix(c(0.1, 0.0, 0.0, 0.0,
                            0.5, 0.4, 0.0, 0.0,
                            0.0, 0.4, 0.7, 0.0,
                            0.0, 0.0, 0.25, 0.4),
                          nrow = 4, ncol = 4, byrow = TRUE)
  
  reproduction_matrix <- matrix(c(0.0, 0.0, 0.0, 2.6,
                                  0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0),
                                nrow = 4, ncol = 4, byrow = TRUE)
  
  # Define clutch distribution
  clutch_distribution <- c(0.06, 0.08, 0.13, 0.15, 0.16, 0.18, 0.15, 0.06, 0.03)
  
  # Define simulation parameters
  num_stages <- 4
  simulation_length <- 24
  initial_population_adult <- c(0, 0, 0, 100)  # All individuals start as adults
  initial_population_spread <- state_initialise_spread(num_stages, 100)  # Spread individuals across life stages
  
  # Run stochastic simulations
  pop_adult <- stochastic_simulation(initial_population_adult, growth_matrix, reproduction_matrix, clutch_distribution, simulation_length)
  pop_spread <- stochastic_simulation(initial_population_spread, growth_matrix, reproduction_matrix, clutch_distribution, simulation_length)
  
  # Create and save the stochastic simulation plot
  png(filename="question_2.png", width = 600, height = 400)
  plot(0:simulation_length, pop_adult, type = "l", col = "red", lwd = 2, ylim = range(c(pop_adult, pop_spread)),
       xlab = "Time (Years)", ylab = "Total Population Size", main = "Stochastic Population Growth Comparison")
  lines(0:simulation_length, pop_spread, col = "blue", lwd = 2)
  legend("topright", legend = c("All Adults", "Evenly Spread"), col = c("red", "blue"), lwd = 2)
  Sys.sleep(0.1)  # Ensure file is saved before closing
  dev.off()
  
  # Explanation of results
  return("In deterministic simulations, the population size changes smoothly over time due to fixed transition rates. However, in stochastic simulations, randomness is introduced at each step, leading to fluctuations in population size. This variation occurs because individual births, deaths, and transitions are probabilistic rather than fixed, making the population growth less predictable and more variable.")
}

# Questions 3 and 4 involve writing code elsewhere to run your simulations on the cluster

# Question 5
question_5 <- function() {
  library(ggplot2)
  
  # Define the four initial condition categories
  init_conditions <- c("big_adult", "small_adult", "big_spread", "small_spread")
  
  # Initialize extinction counters for each category
  extinction_counts <- setNames(rep(0, length(init_conditions)), init_conditions)
  total_counts <- setNames(rep(0, length(init_conditions)), init_conditions)
  
  # Define the directory containing the .rda files
  rda_folder <- "../results"
  
  #  Loop through 100 simulation result files
  for (iter in 1:100) {
    file_name <- paste0(rda_folder, "/demographic_cluster_", iter, ".rda")
    
    # Check if the file exists before attempting to load it
    if (file.exists(file_name)) {
      load(file_name)  # Load the .rda file (should contain `simulation_results`)
      
      # Ensure `simulation_results` exists to prevent errors
      if (!exists("simulation_results")) next
      
      # Correctly map `initial_condition` based on `iter` value
      condition <- ifelse(iter < 26, "big_adult",
                          ifelse(iter < 51, "small_adult",
                                 ifelse(iter < 76, "big_spread", "small_spread")))
      
      # Debugging: Check if condition assignment is correct
      print(paste("iter:", iter, "-> condition:", condition))
      
      # Count total simulations for this condition
      total_counts[condition] <- total_counts[condition] + length(simulation_results)
      
      # Count extinctions (final population size == 0)
      extinction_counts[condition] <- extinction_counts[condition] + sum(sapply(simulation_results, function(sim) {
        if (is.numeric(sim)) return(tail(sim, 1) == 0) else return(NA)
      }), na.rm = TRUE)  # Ignore NA values
    }
  }
  
  #  Compute extinction proportions (avoid division by zero)
  extinction_proportions <- extinction_counts / ifelse(total_counts == 0, 1, total_counts)
  
  # Fix `initial_condition` to ensure correct classification
  df <- data.frame(
    Initial_Condition = factor(names(extinction_proportions), levels = init_conditions),  # Ensures correct order
    Extinction_Probability = extinction_proportions
  )
  
  # Open PNG device for saving the plot
  png("../results/question_5.png", width = 600, height = 400, units = "px", bg = "white")
  
  # nsure X-axis labels display correctly
  extinction_plot <- ggplot(df, aes(x = Initial_Condition, y = Extinction_Probability, fill = Initial_Condition)) +
    geom_bar(stat = "identity", color = "black") +
    theme_minimal(base_size = 14) +
    labs(title = "Extinction Probability by Initial Condition",
         x = "Initial Condition",
         y = "Proportion of Simulations Resulting in Extinction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray95"))
  
  # Save the plot to file
  print(extinction_plot)
  dev.off()  # Close the PNG device
  
  # Identify the population with the highest extinction probability
  most_extinct <- names(which.max(extinction_proportions))
  
  #  Return an explanation
  explanation <- sprintf(
    "The most extinction-prone population was %s. This outcome is expected, as smaller populations are more vulnerable to stochastic fluctuations. 
    In contrast, larger or well-distributed populations tend to be more stable due to the buffering effects provided by greater numbers of individuals and life stages.",
    most_extinct
  ) 
  

  
  return(explanation)
}

# Question 6


# Question 7
species_richness <- function(community) {

  unique_species <- unique(community)
  richness <- length(unique_species) #SR refers to  the number of different species in a sample or community
  return(richness)
}


# Question 8
init_community_max <- function(size){
  community <- seq(1, size)
  return(community) 
}


# Question 9

init_community_min <- function(size) {
  community <- rep(1, size)
  return(community)
}



# Question 10
choose_two <- function(max_value) {
  choose_numbers <- sample(1:max_value, size = 2, replace = FALSE)
  return(choose_numbers)
}
  
  
# Question 11
neutral_step <- function(community) {
  indices <- choose_two(length(community))
  community[indices[1]] <- community[indices[2]]
  return(community)
  }


# Question 12
neutral_generation <- function(community) {
  num_individuals <- length(community)
  num_steps <- round(num_individuals / 2 + sample(c(-0.5, 0.5), 1))
  for (i in 1:num_steps) {
    community <- neutral_step(community)
  }

  return(community)
}


# Question 13
neutral_time_series <- function(community,duration)  {
  # Initialize the species richness time series vector
  richness_series <- numeric(duration + 1)  
  richness_series[1] <- length(unique(community))  # Initial species richness
  
  # Loop through each generation
  for (gen in 1:duration) {
    community <- neutral_generation(community)  # Simulate one neutral generation
    richness_series[gen + 1] <- length(unique(community))  # Record species richness
  }
  
  return(richness_series)
}

# Question 14
question_14 <- function() {
  # Define initial community with maximal diversity
  community <- init_community_max(100)
  
  # Run neutral theory simulation for 200 generations
  richness_series <- neutral_time_series(community, duration = 200)
  
  # Ensure correct filename with .png extension
  png(filename="question_14.png", width = 600, height = 400)
  
  # Plot the species richness time series
  plot(richness_series, type = "l", col = "blue", lwd = 2,
       xlab = "Generation", ylab = "Species Richness",
       main = "Neutral Model: Species Richness Over Time")
  
  Sys.sleep(0.1)  # Allow time for the file to save correctly
  dev.off()  # Close the graphics device to ensure proper saving
  
  return("The system will always converge to a single surviving species (species richness = 1) if we wait long enough. 
          This occurs because, under neutral theory, all species have an equal probability of extinction due to stochastic drift. 
          Eventually, one species will dominate purely by chance, leading to species fixation.")
}


# Question 15
neutral_step_speciation <- function(community, speciation_rate) {
  # Ensure speciation_rate is within a valid range
  if (speciation_rate < 0 || speciation_rate > 1) {
    stop("speciation_rate must be between 0 and 1")
  }
  
  # Select an individual to be replaced
  death_index <- sample(seq_along(community), 1) 
  
  # Speciation occurs with probability speciation_rate
  if (runif(1) < speciation_rate) {
    # Assign a new species number, which must be unique
    existing_species <- unique(community)
    new_species <- max(existing_species) + 1
    
    # Replace the deceased individual with a new species
    community[death_index] <- new_species
  } else {
    # Replace the individual according to the neutral model
    birth_index <- sample(seq_along(community), 1) # Select an individual to reproduce
    community[death_index] <- community[birth_index]
  }
  
  return(community)
}


# Question 16
neutral_generation_speciation <- function(community, speciation_rate) {
  # Ensure speciation_rate is within a valid range
  if (speciation_rate < 0 || speciation_rate > 1) {
    stop("speciation_rate must be between 0 and 1")
  }
  
  # One generation consists of N updates (N is the community size)
  N <- length(community)
  
  # Perform N neutral steps, each with speciation
  for (i in 1:N) {
    community <- neutral_step_speciation(community, speciation_rate)
  }
  
  return(community)
}


# Question 17
neutral_time_series_speciation <- function(initial_community, speciation_rate, duration) {
  # Ensure speciation_rate is within a valid range
  if (speciation_rate < 0 || speciation_rate > 1) {
    stop("speciation_rate must be between 0 and 1")
  }
  
  # Initialize time series with the species richness of the initial community
  species_richness_series <- numeric(duration + 1)
  species_richness_series[1] <- length(unique(initial_community))
  
  # Set the initial community
  community <- initial_community
  
  # Simulate over generations
  for (gen in 1:duration) {
    community <- neutral_generation_speciation(community, speciation_rate)
    species_richness_series[gen + 1] <- length(unique(community))
  }
  
  return(species_richness_series)
}

# Question 18 
question_18 <- function() {
  # Set parameters
  speciation_rate <- 0.1
  community_size <- 100
  duration <- 200
  
  # Initialize two different starting communities
  init_community_max <- seq(1, community_size)  # Each individual is a different species
  init_community_min <- rep(1, community_size)  # All individuals belong to the same species
  
  # Run the simulation
  time_series_max <- neutral_time_series_speciation(init_community_max, speciation_rate, duration)
  time_series_min <- neutral_time_series_speciation(init_community_min, speciation_rate, duration)
  
  # Generate and save the plot
  png(filename = "question_18.png", width = 600, height = 400)
  plot(
    0:duration, time_series_max, type = "l", col = "blue", lwd = 2, ylim = c(0, community_size),
    xlab = "Generations", ylab = "Species Richness",
    main = "Species Richness Over Time with Speciation"
  )
  lines(0:duration, time_series_min, col = "red", lwd = 2)
  legend("topright", legend = c("Max Initial Richness", "Min Initial Richness"), col = c("blue", "red"), lty = 1, lwd = 2)
  
  Sys.sleep(0.1)  # Ensure the plot is fully saved
  dev.off()
  
  # Return explanation
  return("The plot shows that regardless of the initial species richness, the neutral model converges to a similar equilibrium species richness over time. 
  The initially species-rich community (blue) and the initially species-poor community (red) both stabilize at the same species richness after enough generations. 
  This happens because the neutral model assumes all species are ecologically equivalent, leading to a balance between speciation and extinction over time. 
  Thus, the long-term species richness is determined primarily by speciation rate and community size, rather than the starting conditions.")
}

# Question 19
species_abundance <- function(community) {
  # Count the number of individuals per species
  abundance_table <- table(community)
  
  # Convert to a sorted vector in descending order
  abundance_vector <- sort(as.numeric(abundance_table), decreasing = TRUE)
  
  return(abundance_vector)
}



# Question 20
octaves <- function(abundance_vector) {
  # Handle empty input
  if (length(abundance_vector) == 0) {
    return(integer(0))
  }
  
  # Compute octave class index
  octave_classes <- floor(log2(abundance_vector)) + 1
  
  # Determine the maximum required octave class
  max_octave <- max(octave_classes)
  
  # Count species in each octave class
  octave_distribution <- tabulate(factor(octave_classes, levels = 1:max_octave))
  
  return(octave_distribution)
}



# Question 21
sum_vect <- function(x, y) {
  # Get the lengths of both vectors
  len_x <- length(x)
  len_y <- length(y)
  
  # Find the maximum length
  max_len <- max(len_x, len_y)
  
  # Pad x with zeros if it is shorter
  if (len_x < max_len) {
    x <- c(x, rep(0, max_len - len_x))
  }
  # Pad y with zeros if it is shorter
  if (len_y < max_len) {
    y <- c(y, rep(0, max_len - len_y))
  }
  
  # Element-wise addition
  return(x + y)
}

# Question 22
question_22 <- function() {
  # Set parameters
  speciation_rate <- 0.1
  community_size <- 100
  burn_in <- 200
  duration <- 2000
  sampling_interval <- 20
  num_samples <- duration / sampling_interval  # Number of sampling points
  
  # Initialize two different communities
  init_community_max <- seq(1, community_size)  # Each individual is a different species
  init_community_min <- rep(1, community_size)  # All individuals belong to the same species
  
  # Burn-in phase
  community_max <- neutral_generation_speciation(init_community_max, speciation_rate, burn_in)
  community_min <- neutral_generation_speciation(init_community_min, speciation_rate, burn_in)
  
  # Store octave distributions
  octaves_max <- list()
  octaves_min <- list()
  
  # Continue simulation for 2000 generations, sampling every 20 generations
  for (i in 1:num_samples) {
    # Run 20 generations
    community_max <- neutral_generation_speciation(community_max, speciation_rate, sampling_interval)
    community_min <- neutral_generation_speciation(community_min, speciation_rate, sampling_interval)
    
    # Compute octave class distribution
    octaves_max[[i]] <- octaves(species_abundance(community_max))
    octaves_min[[i]] <- octaves(species_abundance(community_min))
  }
  
  # Compute mean species abundance in octave classes
  mean_octaves_max <- Reduce(sum_vect, octaves_max) / num_samples
  mean_octaves_min <- Reduce(sum_vect, octaves_min) / num_samples
  
  # Generate and save bar plot
  png(filename = "question_22.png", width = 600, height = 400)
  par(mfrow = c(1, 2))  # Two subplots
  
  barplot(mean_octaves_max, col = "blue", xlab = "Octave Class", ylab = "Mean Species Count",
          main = "Max Initial Richness", ylim = c(0, max(mean_octaves_max, mean_octaves_min)))
  
  barplot(mean_octaves_min, col = "red", xlab = "Octave Class", ylab = "Mean Species Count",
          main = "Min Initial Richness", ylim = c(0, max(mean_octaves_max, mean_octaves_min)))
  
  Sys.sleep(0.1)  # Ensure the file is written completely
  dev.off()
  
  # Return written explanation
  return("The initial condition of the system does not matter in the long run. 
  Both initial conditions converge to similar species abundance distributions after enough generations. 
  This is because the neutral model assumes that all species have the same ecological fitness. 
  As a result, the system gradually reaches a steady state where species abundance is determined primarily by stochastic birth, death, and speciation processes, rather than by the initial species richness.")
}

# Question 23
neutral_cluster_run <- function(speciation_rate, size, wall_time, interval_rich, interval_oct, burn_in_generations, output_file_name) {
  # Function implementation for running cluster-based simulations (to be completed)
}

# Questions 24 and 25 involve writing code elsewhere to run simulations on the cluster

# Question 26 
process_neutral_cluster_results <- function() {
  
  combined_results <- list() # Create a list output to return
  # Save results to an .rda file (to be implemented)
  
}

plot_neutral_cluster_results <- function() {
  # Load combined_results from the .rda file
  
  png(filename="plot_neutral_cluster_results.png", width = 600, height = 400)
  # Plot the results here
  Sys.sleep(0.1)
  dev.off()
  
  return(combined_results)
}

# Challenge questions - these are substantially harder and worth fewer marks.
# Attempt them only if you've completed all the main questions.

# Challenge question A
Challenge_A <- function() {
  
  png(filename="Challenge_A.png", width = 600, height = 400)
  # Plot your graph here
  Sys.sleep(0.1)
  dev.off()
  
}

# Challenge question B
Challenge_B <- function() {
  
  png(filename="Challenge_B.png", width = 600, height = 400)
  # Plot your graph here
  Sys.sleep(0.1)
  dev.off()
  
}

# Challenge question C
Challenge_C <- function() {
  
  png(filename="Challenge_C.png", width = 600, height = 400)
  # Plot your graph here
  Sys.sleep(0.1)
  dev.off()
  
}

# Challenge question D
Challenge_D <- function() {
  
  png(filename="Challenge_D.png", width = 600, height = 400)
  # Plot your graph here
  Sys.sleep(0.1)
  dev.off()
}

# Challenge question E
Challenge_E <- function() {
  
  png(filename="Challenge_E.png", width = 600, height = 400)
  # Plot your graph here
  Sys.sleep(0.1)
  dev.off()
  
  return("Type your written answer here")
}

