##########################
# 0. Environment Cleanup & Package Loading
##########################
rm(list = ls())    

library(dplyr)
library(ggplot2)
library(minpack.lm)  # Used for nonlinear nlsLM
library(tidyr)
library(segmented)   # In case segmented regression is needed
library(MuMIn)       # For AICc calculation
library(readr)       # For reading/writing CSV files
library(broom)       # For glance/tidy functions

##########################
# 1. Data Reading and Initial Cleaning
##########################

# Read raw data
data <- read.csv("../data/LogisticGrowthData.csv")

message("Loaded ", length(colnames(data)), " columns.")
message("Unique response variable units: ", paste(unique(data$PopBio_units), collapse=", "))
message("Unique independent variable units: ", paste(unique(data$Time_units), collapse=", "))

# Generate unique ID (in this example, concatenating 4 fields as a string then converting to factor)
unique_IDs <- paste(data$Species, data$Temp, data$Medium, data$Citation, sep="_")

data_unique_ID <- data %>%
  mutate(
    ID = factor(
      unique_IDs,
      levels = unique(unique_IDs),
      labels = seq_along(unique(unique_IDs))
    )
  ) %>%
  dplyr::select(-Species, -Temp, -Medium, -Citation)

message("Number of unique datasets: ", dplyr::n_distinct(data_unique_ID$ID))

# Check for missing values
num_missing <- sum(is.na(data_unique_ID))
if (num_missing > 0) {
  message("Warning: ", num_missing, " missing values found.")
}

# Basic filtering
cleaned_growth_data <- data_unique_ID %>%
  filter(Time >= 0, PopBio > 0) %>%
  mutate(logPopBio = log(PopBio))

write.csv(cleaned_growth_data, "../data/modified_growth_data.csv", row.names = FALSE)
message("Cleaned data saved to '../data/modified_growth_data.csv'")
