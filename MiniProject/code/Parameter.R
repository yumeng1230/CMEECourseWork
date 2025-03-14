##############################
# extract_model_parameters.R
# Extract parameters from fitted models
##############################

# 1. Environment Setup & Load Packages
rm(list = ls())
library(dplyr)
library(ggplot2)
library(minpack.lm)  # For nlsLM
library(tidyr)
library(segmented)
library(MuMIn)       # For AICc calculation
library(readr)       # For reading/writing CSV files
library(broom)       # For tidy/glance functions

# 2. Read the cleaned data
data <- read_csv("../data/modified_growth_data.csv")

# 3. Define model functions

# 3.1 Linear Models (Polynomial)
fit_linear_models <- function(df) {
  # Quadratic polynomial
  fit_quad <- lm(logPopBio ~ Time + I(Time^2), data = df)
  # Cubic polynomial
  fit_cubic <- lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df)
  
  list(Quadratic = fit_quad, Cubic = fit_cubic)
}

# 3.2 Nonlinear Model Function Definitions
logistic_model <- function(t, r_max, K, N_0) {
  log(N_0 * K * exp(r_max * t) / (K + N_0 * (exp(r_max * t) - 1)))
}

gompertz_model <- function(t, r_max, K, N_0, t_lag) {
  N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t) / ((K - N_0) * log(10)) + 1))
}

baranyi_model <- function(t, r_max, K, N_0, t_lag) {
  v <- r_max
  m <- 1
  h_0 <- t_lag * r_max
  B <- t + (1 / r_max) * log(exp(-v * t) + exp(-h_0) - exp((-v * t) - h_0))
  N_0 + r_max * B - (1 / m) * log(1 + ((exp(m * r_max * B) - 1) / exp(m * (K - N_0))))
}

three_phase_linear_model <- function(t, N0, Nmax, tLAG, tMAX) {
  ifelse(t <= tLAG,
         N0,
         ifelse(t >= tMAX,
                Nmax,
                N0 + (Nmax - N0) / (tMAX - tLAG) * (t - tLAG)))
}

# 3.3 Multi-Start Random Initialization for nlsLM
compute_max_slope <- function(df) {
  df_sorted <- df %>% arrange(Time)
  slopes <- diff(df_sorted$logPopBio) / diff(df_sorted$Time)
  slopes <- slopes[!is.infinite(slopes) & !is.na(slopes)]
  if (length(slopes) == 0) return(0.1)
  max(slopes)
}

set.seed(123)  # Fix random seed
generate_random_inits <- function(n_starts, df) {
  smax <- compute_max_slope(df)
  if (smax < 0.01) smax <- 0.01
  rmax_min <- 0.01
  rmax_max <- min(smax, 5)
  
  # Adjust K range by scaling down/up the observed logPopBio values
  K_min <- min(df$logPopBio, na.rm = TRUE) * 0.8
  K_max <- max(df$logPopBio, na.rm = TRUE) * 1.2
  
  N0_min <- min(df$logPopBio, na.rm = TRUE)
  N0_max <- max(df$logPopBio, na.rm = TRUE)
  
  tlag_min <- 0
  tlag_max <- max(df$Time, na.rm = TRUE)
  
  tibble(
    r_max = runif(n_starts, rmax_min * 1.01, rmax_max * 0.99),
    K     = runif(n_starts, K_min * 1.01,   K_max * 0.99),
    N_0   = runif(n_starts, N0_min * 1.01,  N0_max * 0.99),
    t_lag = runif(n_starts, tlag_min * 1.01, tlag_max * 0.99)
  )
}

# 3.4 Multi-Start Fitting for Nonlinear Models
multi_start_fit_one_model <- function(df, model_name, n_starts = 150) {
  formula_list <- list(
    "Logistic" = logPopBio ~ logistic_model(Time, r_max, K, N_0),
    "Gompertz" = logPopBio ~ gompertz_model(Time, r_max, K, N_0, t_lag),
    "Baranyi"  = logPopBio ~ baranyi_model(Time, r_max, K, N_0, t_lag)
  )
  fit_formula <- formula_list[[model_name]]
  
  param_names <- if (model_name == "Logistic") {
    c("r_max", "K", "N_0")
  } else {
    c("r_max", "K", "N_0", "t_lag")
  }
  
  inits_tbl <- generate_random_inits(n_starts, df)
  results_list <- list()
  
  for (i in seq_len(n_starts)) {
    start_i <- inits_tbl[i, param_names, drop = FALSE]
    fit_i <- tryCatch(
      nlsLM(fit_formula, data = df, start = as.list(start_i)),
      error = function(e) NULL
    )
    if (!is.null(fit_i)) {
      results_list[[length(results_list) + 1]] <- fit_i
    }
  }
  if (length(results_list) == 0) return(NULL)
  
  best_fit <- NULL
  best_aicc <- Inf
  for (mod in results_list) {
    current_aicc <- AICc(mod)
    if (current_aicc < best_aicc) {
      best_aicc <- current_aicc
      best_fit <- mod
    }
  }
  return(best_fit)
}

# Multi-Start Fitting for the Three-Phase Linear Model
multi_start_fit_three_phase_model <- function(df, n_starts = 150) {
  three_phase_formula <- logPopBio ~ three_phase_linear_model(Time, N0, Nmax, tLAG, tMAX)
  
  t_min <- min(df$Time)
  t_max <- max(df$Time)
  log_min <- min(df$logPopBio)
  log_max <- max(df$logPopBio)
  
  inits_list <- replicate(n_starts, {
    t_rand <- sort(runif(2, min = t_min, max = t_max))
    n_rand <- sort(runif(2, min = log_min, max = log_max))
    c(N0 = n_rand[1], Nmax = n_rand[2], tLAG = t_rand[1], tMAX = t_rand[2])
  })
  inits_mat <- t(inits_list)
  
  results_list <- list()
  for (i in seq_len(n_starts)) {
    start_i <- as.list(inits_mat[i, ])
    fit_i <- tryCatch({
      nlsLM(three_phase_formula, data = df, start = start_i)
    }, error = function(e) NULL)
    if (!is.null(fit_i)) {
      results_list[[length(results_list) + 1]] <- fit_i
    }
  }
  if (length(results_list) == 0) return(NULL)
  
  best_fit <- NULL
  best_aicc <- Inf
  for (mod in results_list) {
    current_aicc <- AICc(mod)
    if (current_aicc < best_aicc) {
      best_aicc <- current_aicc
      best_fit <- mod
    }
  }
  return(best_fit)
}

# Comprehensive fitting of all nonlinear models
fit_nonlinear_models <- function(df, n_starts = 150) {
  if (nrow(df) < 5) return(NULL)
  
  logistic_fit  <- multi_start_fit_one_model(df, "Logistic",  n_starts)
  gompertz_fit  <- multi_start_fit_one_model(df, "Gompertz",  n_starts)
  baranyi_fit   <- multi_start_fit_one_model(df, "Baranyi",   n_starts)
  three_phase_fit <- multi_start_fit_three_phase_model(df, n_starts)
  
  list(
    Logistic   = logistic_fit,
    Gompertz   = gompertz_fit,
    Baranyi    = baranyi_fit,
    ThreePhase = three_phase_fit
  )
}

# 4. Extract model parameters for each ID
linear_params_list <- list()
nonlinear_params_list <- list()

ids <- unique(data$ID)
for (id in ids) {
  df_sub <- data %>% filter(ID == id)
  if(nrow(df_sub) < 5) next
  
  # 4.1 Extract parameters from linear models
  linear_fits <- fit_linear_models(df_sub)
  for(model_name in names(linear_fits)) {
    fit_obj <- linear_fits[[model_name]]
    params <- tidy(fit_obj)
    params <- params %>% mutate(ID = id, Model = model_name)
    linear_params_list[[length(linear_params_list) + 1]] <- params
  }
  
  # 4.2 Extract parameters from nonlinear models
  nonlinear_fits <- fit_nonlinear_models(df_sub, n_starts = 150)
  if (is.null(nonlinear_fits)) next
  for(model_name in names(nonlinear_fits)) {
    fit_obj <- nonlinear_fits[[model_name]]
    if (!is.null(fit_obj)) {
      params <- tidy(fit_obj)
      params <- params %>% mutate(ID = id, Model = model_name)
      nonlinear_params_list[[length(nonlinear_params_list) + 1]] <- params
    }
  }
}

# Combine the parameters from both linear and nonlinear models
linear_params_df <- bind_rows(linear_params_list)
nonlinear_params_df <- bind_rows(nonlinear_params_list)
all_params_df <- bind_rows(linear_params_df, nonlinear_params_df)

# 5. Save the results to a CSV file
write_csv(all_params_df, "../results/model_parameters_extracted.csv")
message("Extracted model parameters have been saved to '../results/model_parameters_extracted.csv'.")
