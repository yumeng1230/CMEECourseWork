
##########################
# 2. Linear (Polynomial) Regression Comparison
##########################

library(dplyr)
library(broom)
library(MuMIn)  # For calculating AICc

# Read data
datalinear <- read.csv("../data/modified_growth_data.csv")

# Function: Fit quadratic & cubic polynomial on a data subset and output the results
fit_models <- function(df) {
  if (nrow(df) < 8) return(NULL)  # Filter out datasets that are too small
  
  # Quadratic polynomial
  model1 <- lm(logPopBio ~ Time + I(Time^2), data = df)
  summary1 <- glance(model1)
  coef1 <- tidy(model1) %>% filter(term == "Time") %>% pull(estimate)
  
  # Cubic polynomial
  model2 <- lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df)
  summary2 <- glance(model2)
  coef2 <- tidy(model2) %>% filter(term == "Time") %>% pull(estimate)
  
  tibble(
    ID = unique(df$ID),
    Model1_AIC = AIC(model1),
    Model1_AICc = AICc(model1),
    Model1_BIC = BIC(model1),
    Model1_R2 = summary1$r.squared,
    Model1_B1 = coef1,
    
    Model2_AIC = AIC(model2),
    Model2_AICc = AICc(model2),
    Model2_BIC = BIC(model2),
    Model2_R2 = summary2$r.squared,
    Model2_B1 = coef2,
    
    Best_Model = ifelse(AICc(model1) < AICc(model2), "Quadratic", "Cubic")
  )
}

# Group by each ID and fit models
results_linear <- datalinear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_models) %>%
  bind_rows()

# Select only columns with AICc and BIC
results_selected <- results_linear %>%
  dplyr::select(ID, Model1_AICc, Model1_BIC, Model2_AICc, Model2_BIC)

# Save full results
write.csv(results_linear, "../results/model_comparison.csv", row.names = FALSE)
message("Full model comparison results saved to '../results/model_comparison.csv'")

# Save only AICc and BIC results
write.csv(results_selected, "../results/model_aicc_bic_comparison.csv", row.names = FALSE)
message("AICc and BIC comparison saved to '../results/model_aicc_bic_comparison.csv'")

##########################
# 3. Define Nonlinear Growth Models
##########################

logistic_model <- function(t, r_max, K, N_0) {
  log(N_0 * K * exp(r_max * t) / (K + N_0 * (exp(r_max * t) - 1)))
}

gompertz_model <- function(t, r_max, K, N_0, t_lag) {
  N_0 + (K - N_0) * exp(
    -exp(
      r_max * exp(1) * (t_lag - t) / ((K - N_0) * log(10)) + 1
    )
  )
}

baranyi_model <- function(t, r_max, K, N_0, t_lag) {
  v <- r_max
  m <- 1
  h_0 <- t_lag * r_max
  B <- t + (1 / r_max) * log(exp(-v * t) + exp(-h_0) - exp((-v * t) - h_0))
  
  N_0 + r_max * B - (1 / m) * log(
    1 + ((exp(m * r_max * B) - 1) / exp(m * (K - N_0)))
  )
}

##########################
# 3.1 Define "Three-Phase Linear Model" (New)
##########################
# According to the problem description, the model expressed on a logarithmic scale:
#   For t <= tLAG,         N(t) = N0
#   For tLAG < t < tMAX,    N(t) = N0 + (Nmax - N0)/(tMAX - tLAG) * (t - tLAG)
#   For t >= tMAX,          N(t) = Nmax
three_phase_linear_model <- function(t, N0, Nmax, tLAG, tMAX) {
  ifelse(t <= tLAG,
         N0,
         ifelse(t >= tMAX,
                Nmax,
                N0 + (Nmax - N0)/(tMAX - tLAG) * (t - tLAG)
         )
  )
}

############################## 
# 4. Multi-Start Random Initialization  
##############################
compute_max_slope <- function(df) {
  df_sorted <- df %>% arrange(Time)
  slope_vals <- diff(df_sorted$logPopBio) / diff(df_sorted$Time)
  slope_vals <- slope_vals[!is.infinite(slope_vals) & !is.na(slope_vals)]
  if (length(slope_vals) == 0) return(0.1)
  max(slope_vals)
}

set.seed(123)  # Set random seed

generate_random_inits <- function(n_starts, df) {
  smax <- compute_max_slope(df)
  if (smax < 0.01) smax <- 0.01
  
  rmax_min <- 0.01
  rmax_max <- min(smax, 5)  # Limit r_max to a maximum of 5
  
  K_min <- min(df$logPopBio, na.rm = TRUE) 
  K_max <- max(df$logPopBio, na.rm = TRUE)   # Provide a more appropriate range for K
  
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

##########################
# 5. Batch Nonlinear Fitting & Saving
##########################

# Read the cleaned data
data_non_linear <- read.csv("../data/modified_growth_data.csv")

# Original function: multi_start_fit_one_model
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
    start_i <- inits_tbl[i, param_names, drop=FALSE]
    fit_i <- tryCatch(
      {
        nlsLM(fit_formula, data = df, start = as.list(start_i))
      },
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

# Original function: fit_nonlinear_models
fit_nonlinear_models <- function(df, n_starts = 150) {
  if (nrow(df) < 8) return(NULL)
  
  logistic_fit  <- multi_start_fit_one_model(df, "Logistic",  n_starts)
  gompertz_fit  <- multi_start_fit_one_model(df, "Gompertz",  n_starts)
  baranyi_fit   <- multi_start_fit_one_model(df, "Baranyi",   n_starts)
  
  get_model_metrics <- function(model) {
    if (is.null(model)) return(c(NA, NA, NA, NA))
    g <- glance(model)
    c(AIC(model), AICc(model), BIC(model), g$r.squared)
  }
  
  logistic_metrics  <- get_model_metrics(logistic_fit)
  gompertz_metrics  <- get_model_metrics(gompertz_fit)
  baranyi_metrics   <- get_model_metrics(baranyi_fit)
  
  all_aicc <- c(logistic_metrics[2], gompertz_metrics[2], baranyi_metrics[2])
  model_names <- c("Logistic", "Gompertz", "Baranyi")
  best_model <- model_names[which.min(all_aicc)]
  
  tibble(
    ID = unique(df$ID),
    Logistic_AIC  = logistic_metrics[1],
    Logistic_AICc = logistic_metrics[2],
    Logistic_BIC  = logistic_metrics[3],
    Logistic_R2   = logistic_metrics[4],
    
    Gompertz_AIC  = gompertz_metrics[1],
    Gompertz_AICc = gompertz_metrics[2],
    Gompertz_BIC  = gompertz_metrics[3],
    Gompertz_R2   = gompertz_metrics[4],
    
    Baranyi_AIC   = baranyi_metrics[1],
    Baranyi_AICc  = baranyi_metrics[2],
    Baranyi_BIC   = baranyi_metrics[3],
    Baranyi_R2    = baranyi_metrics[4],
    
    Best_Model = best_model
  )
}

# Batch fit the original 3 nonlinear models
results_nonlinear <- data_non_linear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_nonlinear_models, n_starts = 150) %>%
  bind_rows()

# Save nonlinear model comparison results
write.csv(results_nonlinear, "../results/nonlinear_model_comparison.csv", row.names = FALSE)
message("Nonlinear model comparison (multi-start) saved.")

##############################
# 5.1 New Function: Multi-Start Fitting for Three-Phase Linear Model
##############################

multi_start_fit_three_phase_model <- function(df, n_starts = 150) {
  # Define the formula specific to the three-phase linear model
  three_phase_formula <- logPopBio ~ three_phase_linear_model(Time, N0, Nmax, tLAG, tMAX)
  
  # Set range for initial values
  df_sorted <- df %>% arrange(Time)
  t_min <- min(df$Time)
  t_max <- max(df$Time)
  
  log_min <- min(df$logPopBio)
  log_max <- max(df$logPopBio)
  
  # Generate random initial values
  inits_list <- replicate(n_starts, {
    t_rand <- sort(runif(2, min = t_min, max = t_max))
    n_rand <- sort(runif(2, min = log_min, max = log_max))
    
    c(N0 = n_rand[1],
      Nmax = n_rand[2],
      tLAG = t_rand[1],
      tMAX = t_rand[2])
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

# New function: Batch fit the three-phase linear model and output a results table
fit_three_phase_linear_models <- function(df, n_starts = 150) {
  if (nrow(df) < 8) return(NULL)
  
  three_phase_fit <- multi_start_fit_three_phase_model(df, n_starts)
  
  get_model_metrics <- function(model) {
    if (is.null(model)) return(c(NA, NA, NA, NA))
    g <- glance(model)
    c(AIC(model), AICc(model), BIC(model), g$r.squared)
  }
  three_phase_metrics <- get_model_metrics(three_phase_fit)
  
  tibble(
    ID = unique(df$ID),
    ThreePhase_AIC  = three_phase_metrics[1],
    ThreePhase_AICc = three_phase_metrics[2],
    ThreePhase_BIC  = three_phase_metrics[3],
    ThreePhase_R2   = three_phase_metrics[4]
  )
}

# Batch fit the three-phase linear model for all IDs
results_three_phase <- data_non_linear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_three_phase_linear_models, n_starts = 150) %>%
  bind_rows()

# Save three-phase linear model comparison results
write.csv(results_three_phase, "../results/three_phase_model_comparison.csv", row.names = FALSE)
message("Three-phase linear model comparison saved.")

##########################
# 6. Final Merge & Global Best Model Selection based on AICc
##########################

library(dplyr)
library(readr)

# Read linear model results
linear_results <- read_csv("../results/model_comparison.csv") %>%
  dplyr::select(ID, Model1_AICc, Model1_BIC, Model2_AICc, Model2_BIC) %>%
  rename(
    Quadratic_AICc = Model1_AICc,
    Cubic_AICc     = Model2_AICc,
    Quadratic_BIC  = Model1_BIC,
    Cubic_BIC      = Model2_BIC
  )

# Read nonlinear model results (Logistic/Gompertz/Baranyi)
nonlinear_results <- read_csv("../results/nonlinear_model_comparison.csv") %>%
  dplyr::select(
    ID,
    Logistic_AICc, Gompertz_AICc, Baranyi_AICc,
    Logistic_BIC,  Gompertz_BIC,  Baranyi_BIC
  )

# Read three-phase linear model results
three_phase_results <- read_csv("../results/three_phase_model_comparison.csv") %>%
  dplyr::select(
    ID,
    ThreePhase_AICc, ThreePhase_BIC
  )

# Merge linear, original 3 nonlinear, and new three-phase model results
final_results <- linear_results %>%
  full_join(nonlinear_results, by = "ID") %>%
  full_join(three_phase_results, by = "ID")

write_csv(final_results, "../results/final_model_comparison.csv")

# Read the merged results
df <- read.csv("../results/final_model_comparison.csv", stringsAsFactors = FALSE)

# Replace NA with a large value (1e6) to avoid NA being selected as the best model
df[is.na(df)] <- 1e6

# Define columns for AICc and BIC (now there are 6 models)
model_aicc_cols <- c(
  "Quadratic_AICc", "Cubic_AICc",
  "Logistic_AICc", "Gompertz_AICc", "Baranyi_AICc",
  "ThreePhase_AICc"
)
model_bic_cols  <- c(
  "Quadratic_BIC",  "Cubic_BIC",
  "Logistic_BIC",   "Gompertz_BIC", "Baranyi_BIC",
  "ThreePhase_BIC"
)

# Select the model with the smallest AICc
best_model_aicc <- apply(df[, model_aicc_cols], 1, function(x) {
  c("Quadratic", "Cubic", "Logistic", "Gompertz", "Baranyi", "ThreePhase")[which.min(x)]
})

# Select the model with the smallest BIC
best_model_bic <- apply(df[, model_bic_cols], 1, function(x) {
  c("Quadratic", "Cubic", "Logistic", "Gompertz", "Baranyi", "ThreePhase")[which.min(x)]
})

# Add best model columns
df$Best_Model_AICc <- best_model_aicc
df$Best_Model_BIC  <- best_model_bic

# Save final best model selection results
write.csv(df, "../results/best_model_comparison_updated.csv", row.names = FALSE)
message("Final best model selection (AICc & BIC) saved to '../results/best_model_comparison_updated.csv'.")

# Print first few rows for inspection
cat("===== HEAD OF df =====\n")
print(head(df))

# Count best model selections (AICc)
best_model_counts_aicc <- table(df$Best_Model_AICc)
message("Best Model Selection Counts (AICc):")
print(best_model_counts_aicc)

# Count best model selections (BIC)
best_model_counts_bic <- table(df$Best_Model_BIC)
message("Best Model Selection Counts (BIC):")
print(best_model_counts_bic)

# Count the number of times each model's AICc is set to 1e6
model_onee6_counts_aicc <- sapply(model_aicc_cols, function(col) sum(df[[col]] == 1e6))
message("Each model's AICc is 1e6 for these many rows:")
print(model_onee6_counts_aicc)

# Count the number of times each model's BIC is set to 1e6
model_onee6_counts_bic <- sapply(model_bic_cols, function(col) sum(df[[col]] == 1e6))
message("Each model's BIC is 1e6 for these many rows:")
print(model_onee6_counts_bic)

################################
# Plot Bar Chart: Compare AICc & BIC
################################
library(ggplot2)
library(dplyr)
library(readr)

# Read final best model selection results
model_selection_data <- read_csv("../results/best_model_comparison_updated.csv")

# Count the number of times each model is selected by AICc and BIC
aicc_counts <- table(model_selection_data$Best_Model_AICc)
bic_counts  <- table(model_selection_data$Best_Model_BIC)

# Convert to data frame for ggplot
model_counts_df <- data.frame(
  Model     = rep(names(aicc_counts), 2),
  Count     = c(as.numeric(aicc_counts), as.numeric(bic_counts)),
  Criterion = rep(c("AICc", "BIC"), each = length(aicc_counts))
)

# Create bar plot object
p <- ggplot(model_counts_df, aes(x = Model, y = Count, fill = Criterion)) +
  geom_bar(stat = "identity", position = "dodge") +  # Bars displayed side-by-side
  labs(title = "Comparison of Model Selection by AICc and BIC",
       x = "Models", y = "Selection Counts") +
  theme_minimal() +
  # The color scale is set here; remove scale_fill_manual(values=...) if you prefer default colors
  scale_fill_manual(values = c("AICc" = "#1f77b4", "BIC" = "#ff7f0e")) +
  geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.9))

# Save the plot
ggsave("../results/model_selection_comparison.png", plot = p, width = 8, height = 6, dpi = 300)
message("Plot saved successfully to '../results/model_selection_comparison.png'")



############ Below is an additional example: Plot the fitted curves of 6 models for each ID ############
# This section extends the previous "5 models" plotting to 6 models
# Following the "do not modify original code" principle, we created a new function fit_and_select_best_model_6()
#####################################################################

# Read data
growth_data <- read_csv("../data/modified_growth_data.csv")
all_ids <- sort(unique(growth_data$ID))

# The original fit_and_select_best_model_5() remains unchanged; here we write a new 6-model version:
fit_and_select_best_model_6 <- function(df_sub, n_starts = 50) {
  
  # 1) Quadratic
  fit_quad <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2), data = df_sub),
    error = function(e) NULL
  )
  # 2) Cubic
  fit_cubic <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df_sub),
    error = function(e) NULL
  )
  # 3) Logistic
  fit_logistic <- tryCatch(
    multi_start_fit_one_model(df_sub, "Logistic", n_starts),
    error = function(e) NULL
  )
  # 4) Gompertz
  fit_gompertz <- tryCatch(
    multi_start_fit_one_model(df_sub, "Gompertz", n_starts),
    error = function(e) NULL
  )
  # 5) Baranyi
  fit_baranyi <- tryCatch(
    multi_start_fit_one_model(df_sub, "Baranyi", n_starts),
    error = function(e) NULL
  )
  # 6) Three-phase linear
  fit_threephase <- tryCatch(
    multi_start_fit_three_phase_model(df_sub, n_starts),
    error = function(e) NULL
  )
  
  # If any model is NULL, it indicates a fitting failure -> return NULL
  if(is.null(fit_quad)      ||
     is.null(fit_cubic)     ||
     is.null(fit_logistic)  ||
     is.null(fit_gompertz)  ||
     is.null(fit_baranyi)   ||
     is.null(fit_threephase)) {
    return(NULL)
  }
  
  # Calculate AICc and find the best model
  model_list <- list(
    Quadratic  = fit_quad,
    Cubic      = fit_cubic,
    Logistic   = fit_logistic,
    Gompertz   = fit_gompertz,
    Baranyi    = fit_baranyi,
    ThreePhase = fit_threephase
  )
  aicc_values <- sapply(model_list, AICc)  # MuMIn::AICc
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  
  # Generate prediction grid
  t_min <- min(df_sub$Time, na.rm=TRUE)
  t_max <- max(df_sub$Time, na.rm=TRUE)
  new_time <- seq(t_min, t_max, length.out = 200)
  
  # Predict for each model
  # Note: For lm objects, use predict() directly; nls objects work similarly with predict()
  preds_list <- list(
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_quad, newdata = data.frame(Time = new_time)),
      Model = "Quadratic"
    ),
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_cubic, newdata = data.frame(Time = new_time)),
      Model = "Cubic"
    ),
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_logistic, newdata = data.frame(Time = new_time)),
      Model = "Logistic"
    ),
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_gompertz, newdata = data.frame(Time = new_time)),
      Model = "Gompertz"
    ),
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_baranyi, newdata = data.frame(Time = new_time)),
      Model = "Baranyi"
    ),
    data.frame(
      Time = new_time,
      logPopBio = predict(fit_threephase, newdata = data.frame(Time = new_time)),
      Model = "ThreePhase"
    )
  )
  
  preds_all <- do.call(rbind, preds_list)
  
  # Return a list containing "predictions" and "best_model"
  return(list(
    predictions = preds_all,
    best_model  = best_model_name
  ))
}


# Batch generate plots
output_dir <- "../results/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for(this_id in all_ids) {
  # Subset the data
  df_sub <- growth_data %>% filter(ID == this_id)
  
  # Skip if there are fewer than 8 data points
  if(nrow(df_sub) < 8) {
    next
  }
  
  # Fit 6 models & compute best model
  fit_result <- fit_and_select_best_model_6(df_sub, n_starts = 150)
  
  # If any model fitting failed (NULL returned) -> skip
  if(is.null(fit_result)) {
    next
  }
  
  # Otherwise, all models fit successfully; get predictions & best model
  pred_df   <- fit_result$predictions
  best_name <- fit_result$best_model
  
  # Plot: original scatter points + 6 curves
  p <- ggplot() + 
    geom_point(data = df_sub, aes(x = Time, y = logPopBio)) + 
    geom_line(data = pred_df, aes(x = Time, y = logPopBio, color = Model), size = 0.8) + 
    labs(
      title = paste0("ID = ", this_id, ", Best Model (AICc) = ", best_name),
      x = "Time(Hours)",
      y = "log(PopBio)"
    ) + 
    theme_minimal() + 
    scale_color_discrete(name = "Model")+ 
    theme(
      legend.text = element_text(size = 14),   # Increase legend text size
      legend.title = element_text(size = 16, face = "bold")  # Increase and bold legend title
    )
  
  # Save PNG
  png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")  
  png(png_filename, width = 2000, height = 1500, res = 200)  
  print(p)  
  dev.off()  
  
  message("PNG output: ", png_filename, " ; Best model = ", best_name) 
}

message("All tasks completed: Now there are comparison results and visualizations for 6 models including the three-phase linear model.")

# Define function to calculate Akaike weights
calculate_akaike_weights <- function(aic_values) {
  # Calculate ΔAICc
  delta_aic <- aic_values - min(aic_values)
  # Calculate relative likelihoods
  rel_likelihood <- exp(-0.5 * delta_aic)
  # Normalize to get weights
  weights <- rel_likelihood / sum(rel_likelihood)
  return(weights)
}

# Assume the final_results data frame contains AICc columns for 6 models
# Column names: Quadratic_AICc, Cubic_AICc, Logistic_AICc, Gompertz_AICc, Baranyi_AICc, ThreePhase_AICc
model_aicc_cols <- c("Quadratic_AICc", "Cubic_AICc", "Logistic_AICc", 
                     "Gompertz_AICc", "Baranyi_AICc", "ThreePhase_AICc")

# For each row (each ID's data), calculate Akaike weights and store the result in a list
akaike_weights_list <- apply(df[, model_aicc_cols], 1, calculate_akaike_weights)

# If needed, convert the result to a data frame (each row for an ID, each column for a model's Akaike weight)
akaike_weights_df <- as.data.frame(t(akaike_weights_list))
colnames(akaike_weights_df) <- c("Quadratic_w", "Cubic_w", "Logistic_w", 
                                 "Gompertz_w", "Baranyi_w", "ThreePhase_w")

# Merge with final results
final_results_with_weights <- cbind(df, akaike_weights_df)

# View results
head(final_results_with_weights)
# Write the results to a new CSV file
write.csv(final_results_with_weights, "../results/final_model_comparison_with_akaike_weights.csv", row.names = FALSE)
message("Final model comparison results with Akaike weights have been saved to '../results/final_model_comparison_with_akaike_weights.csv'")

# Assume final_results_with_weights contains the Akaike weights for each ID,
# with column names: Quadratic_w, Cubic_w, Logistic_w, Gompertz_w, Baranyi_w, ThreePhase_w

# For each ID, select the model with the largest Akaike weight
final_results_with_weights$Best_Model <- apply(
  final_results_with_weights[, c("Quadratic_w", "Cubic_w", "Logistic_w", 
                                 "Gompertz_w", "Baranyi_w", "ThreePhase_w")],
  1,
  function(x) names(x)[which.max(x)]
)

# View best model for each ID
head(final_results_with_weights[, c("ID", "Best_Model")])

# Count the number of times each model is selected as the best
best_model_counts <- table(final_results_with_weights$Best_Model)
print(best_model_counts)

# Save results to a new file
write.csv(final_results_with_weights, "../results/final_model_comparison_with_weights_and_best.csv", row.names = FALSE)
message("File with Akaike weights and best model selection results has been saved.")

# Assume final_results_with_weights already contains the Akaike weights for each ID,
# with column names: Quadratic_w, Cubic_w, Logistic_w, Gompertz_w, Baranyi_w, ThreePhase_w

# Calculate the average Akaike weight for each model
library(dplyr)
summary_weights <- final_results_with_weights %>%
  summarise(
    Quadratic_mean = mean(Quadratic_w, na.rm = TRUE),
    Cubic_mean     = mean(Cubic_w, na.rm = TRUE),
    Logistic_mean  = mean(Logistic_w, na.rm = TRUE),
    Gompertz_mean  = mean(Gompertz_w, na.rm = TRUE),
    Baranyi_mean   = mean(Baranyi_w, na.rm = TRUE),
    ThreePhase_mean = mean(ThreePhase_w, na.rm = TRUE)
  )
print(summary_weights)

# Convert data to long format for plotting
library(tidyr)
library(ggplot2)
weights_long <- final_results_with_weights %>%
  dplyr::select(ID, Quadratic_w, Cubic_w, Logistic_w, Gompertz_w, Baranyi_w, ThreePhase_w) %>%
  tidyr::pivot_longer(cols = -ID, names_to = "Model", values_to = "Weight")

# Plot boxplot of Akaike weights for each model
library(ggplot2)
library(dplyr)

# Compute median for each model (still useful if needed later)
medians <- weights_long %>%
  group_by(Model) %>%
  summarise(median_weight = median(Weight, na.rm = TRUE))

# Create boxplot with jittered scatter points (without annotations)
p <- ggplot(weights_long, aes(x = Model, y = Weight, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +             # Do not show outliers to keep plot clean
  geom_jitter(width = 0.1, alpha = 0.5, color = "black") +  # Add jittered points to show all data
  scale_fill_brewer(palette = "Set3") +          # Use the Set3 color palette
  labs(title = "Distribution of Akaike Weights by Model",
       x = "Model",
       y = "Akaike Weight") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(face = "bold", hjust = 0.5))

# Save the enhanced boxplot
ggsave("../results/akaike_weights_boxplot_nice.png", plot = p, width = 8, height = 6, dpi = 300)

# Save the boxplot
ggsave("../results/akaike_weights_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)
message("Boxplot saved to '../results/akaike_weights_boxplot.png'")

# Determine the best model for each ID based on Akaike weights
final_results_with_weights$Best_Model_Akaike <- apply(
  final_results_with_weights[, c("Quadratic_w", "Cubic_w", "Logistic_w",
                                 "Gompertz_w", "Baranyi_w", "ThreePhase_w")],
  1,
  function(x) names(x)[which.max(x)]
)
# Count the number of times each model is selected as best
best_model_counts <- table(final_results_with_weights$Best_Model_Akaike)
print(best_model_counts)

# Write the final results (including Akaike weights and best model information) to a new file
write.csv(final_results_with_weights, "../results/final_model_comparison_with_akaike_weights_analysis.csv", row.names = FALSE)
message("New file with Akaike weights and analysis results has been saved.")

