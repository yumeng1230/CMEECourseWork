#############################################
# 01_data_preprocessing.R
# 1) Read raw data
# 2) Perform basic cleaning
# 3) (Optional) Fit polynomial (quadratic/cubic) models
# 4) Save cleaned data to ../data/modified_growth_data.csv
#############################################
rm(list = ls())    

library(dplyr)
library(ggplot2)
library(minpack.lm)  # for nonlinear nlsLM
library(tidyr)
library(segmented)   # if segmented regression is needed
library(MuMIn)       # for AICc calculation
library(readr)       # for reading/writing CSV
library(broom)       # for glance, tidy, etc.

##########################
# 1. Data Reading and Initial Cleaning
##########################

# Read the raw data
data <- read.csv("../data/LogisticGrowthData.csv")

message("Loaded ", length(colnames(data)), " columns.")
message("Unique response variable units: ", paste(unique(data$PopBio_units), collapse=", "))
message("Unique independent variable units: ", paste(unique(data$Time_units), collapse=", "))

# Generate unique ID by concatenating four fields and converting to factor
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

# Basic filtering: keep rows with Time >= 0 and PopBio > 0, and compute log(PopBio)
cleaned_data <- data_unique_ID %>%
  filter(Time >= 0, PopBio > 0) %>%
  mutate(logPopBio = log(PopBio))

write.csv(cleaned_data, "../data/modified_growth_data.csv", row.names = FALSE)
message("Cleaned data saved to '../data/modified_growth_data.csv'")

#### (Optional) Polynomial (quadratic & cubic) regression ####
fit_polynomial_models <- function(df){
  if(nrow(df) < 8) return(NULL)
  # Quadratic
  m_quad  <- lm(logPopBio ~ Time + I(Time^2), data=df)
  # Cubic
  m_cubic <- lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data=df)
  tibble(
    ID = unique(df$ID),
    Quadratic_AICc = AICc(m_quad),
    Quadratic_BIC  = BIC(m_quad),
    Cubic_AICc     = AICc(m_cubic),
    Cubic_BIC      = BIC(m_cubic)
  )
}

# Fit polynomials per ID
results_linear <- cleaned_data %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_polynomial_models) %>%
  bind_rows()

# Save polynomial model comparison
write_csv(results_linear, "../results/model_comparison_linear.csv")
message("Polynomial model comparison (Quadratic vs. Cubic) saved: ../results/model_comparison_linear.csv")

