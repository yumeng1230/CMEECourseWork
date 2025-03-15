##########################
# 0. Environment Cleaning & Package Loading
##########################
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
cleaned_growth_data <- data_unique_ID %>%
  filter(Time >= 0, PopBio > 0) %>%
  mutate(logPopBio = log(PopBio))

write.csv(cleaned_growth_data, "../data/modified_growth_data.csv", row.names = FALSE)
message("Cleaned data saved to '../data/modified_growth_data.csv'")

##########################
# 2. Linear (Polynomial) Regression Comparison
##########################

# Read the modified data
datalinear <- read.csv("../data/modified_growth_data.csv")

# Function to fit quadratic and cubic models on a subset and output results
fit_models <- function(df) {
  if(nrow(df) < 8) return(NULL)  # skip very small datasets
  
  # Quadratic model
  model1 <- lm(logPopBio ~ Time + I(Time^2), data = df)
  summary1 <- glance(model1)
  coef1 <- tidy(model1) %>% filter(term == "Time") %>% pull(estimate)
  
  # Cubic model
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

results_linear <- datalinear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_models) %>%
  bind_rows()

results_selected <- results_linear %>%
  dplyr::select(ID, Model1_AICc, Model1_BIC, Model2_AICc, Model2_BIC)

write.csv(results_linear, "../results/model_comparison.csv", row.names = FALSE)
message("Full model comparison results saved to '../results/model_comparison.csv'")
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
# 3.1 Define Three-Phase Linear Model (New)
##########################
# On a log scale, the model is defined as:
#   For t <= tLAG:       N(t) = N0
#   For tLAG < t < tMAX:  N(t) = N0 + (Nmax - N0)/(tMAX - tLAG) * (t - tLAG)
#   For t >= tMAX:       N(t) = Nmax
three_phase_linear_model <- function(t, N0, Nmax, tLAG, tMAX) {
  ifelse(t <= tLAG,
         N0,
         ifelse(t >= tMAX,
                Nmax,
                N0 + (Nmax - N0)/(tMAX - tLAG) * (t - tLAG)
         )
  )
}

##########################
# 4. Multi-start Random Initial Values
##########################

compute_max_slope <- function(df) {
  df_sorted <- df %>% arrange(Time)
  slope_vals <- diff(df_sorted$logPopBio) / diff(df_sorted$Time)
  slope_vals <- slope_vals[!is.infinite(slope_vals) & !is.na(slope_vals)]
  if(length(slope_vals) == 0) return(0.1)
  max(slope_vals)
}

set.seed(123)

generate_random_inits <- function(n_starts, df) {
  smax <- compute_max_slope(df)
  if(smax < 0.01) smax <- 0.01
  
  rmax_min <- 0.01
  rmax_max <- min(smax, 5)
  
  K_min <- min(df$logPopBio, na.rm = TRUE)
  K_max <- max(df$logPopBio, na.rm = TRUE)
  
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

data_non_linear <- read.csv("../data/modified_growth_data.csv")

multi_start_fit_one_model <- function(df, model_name, n_starts = 150) {
  formula_list <- list(
    "Logistic" = logPopBio ~ logistic_model(Time, r_max, K, N_0),
    "Gompertz" = logPopBio ~ gompertz_model(Time, r_max, K, N_0, t_lag),
    "Baranyi"  = logPopBio ~ baranyi_model(Time, r_max, K, N_0, t_lag)
  )
  
  fit_formula <- formula_list[[model_name]]
  
  param_names <- if(model_name == "Logistic") {
    c("r_max", "K", "N_0")
  } else {
    c("r_max", "K", "N_0", "t_lag")
  }
  
  inits_tbl <- generate_random_inits(n_starts, df)
  
  results_list <- list()
  
  for(i in seq_len(n_starts)) {
    start_i <- inits_tbl[i, param_names, drop=FALSE]
    fit_i <- tryCatch({
      nlsLM(fit_formula, data = df, start = as.list(start_i))
    }, error = function(e) NULL)
    if(!is.null(fit_i)) {
      results_list[[length(results_list) + 1]] <- fit_i
    }
  }
  
  if(length(results_list) == 0) return(NULL)
  
  best_fit <- NULL
  best_aicc <- Inf
  for(mod in results_list) {
    current_aicc <- AICc(mod)
    if(current_aicc < best_aicc) {
      best_aicc <- current_aicc
      best_fit <- mod
    }
  }
  return(best_fit)
}

fit_nonlinear_models <- function(df, n_starts = 150) {
  if(nrow(df) < 8) return(NULL)
  
  logistic_fit  <- multi_start_fit_one_model(df, "Logistic", n_starts)
  gompertz_fit  <- multi_start_fit_one_model(df, "Gompertz", n_starts)
  baranyi_fit   <- multi_start_fit_one_model(df, "Baranyi", n_starts)
  
  get_model_metrics <- function(model) {
    if(is.null(model)) return(c(NA, NA, NA, NA))
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

results_nonlinear <- data_non_linear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_nonlinear_models, n_starts = 150) %>%
  bind_rows()

write.csv(results_nonlinear, "../results/nonlinear_model_comparison.csv", row.names = FALSE)
message("Nonlinear model comparison (multi-start) saved.")

##############################
# 5.1 New Function: Multi-start Fitting for Three-Phase Linear Model
##############################

multi_start_fit_three_phase_model <- function(df, n_starts = 150) {
  three_phase_formula <- logPopBio ~ three_phase_linear_model(Time, N0, Nmax, tLAG, tMAX)
  
  df_sorted <- df %>% arrange(Time)
  t_min <- min(df$Time)
  t_max <- max(df$Time)
  
  log_min <- min(df$logPopBio)
  log_max <- max(df$logPopBio)
  
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
  for(i in seq_len(n_starts)) {
    start_i <- as.list(inits_mat[i, ])
    fit_i <- tryCatch({
      nlsLM(three_phase_formula, data = df, start = start_i)
    }, error = function(e) NULL)
    if(!is.null(fit_i)) {
      results_list[[length(results_list) + 1]] <- fit_i
    }
  }
  
  if(length(results_list) == 0) return(NULL)
  
  best_fit <- NULL
  best_aicc <- Inf
  for(mod in results_list) {
    current_aicc <- AICc(mod)
    if(current_aicc < best_aicc) {
      best_aicc <- current_aicc
      best_fit <- mod
    }
  }
  return(best_fit)
}

fit_three_phase_linear_models <- function(df, n_starts = 150) {
  if(nrow(df) < 8) return(NULL)
  three_phase_fit <- multi_start_fit_three_phase_model(df, n_starts)
  
  get_model_metrics <- function(model) {
    if(is.null(model)) return(c(NA, NA, NA, NA))
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

results_three_phase <- data_non_linear %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(fit_three_phase_linear_models, n_starts = 150) %>%
  bind_rows()

write.csv(results_three_phase, "../results/three_phase_model_comparison.csv", row.names = FALSE)
message("Three-phase linear model comparison saved.")

##########################
# 6. Merge Results & Choose Global Best AICc
##########################

linear_results <- read_csv("../results/model_comparison.csv") %>%
  dplyr::select(ID, Model1_AICc, Model1_BIC, Model2_AICc, Model2_BIC) %>%
  rename(
    Quadratic_AICc = Model1_AICc,
    Cubic_AICc     = Model2_AICc,
    Quadratic_BIC  = Model1_BIC,
    Cubic_BIC      = Model2_BIC
  )

nonlinear_results <- read_csv("../results/nonlinear_model_comparison.csv") %>%
  dplyr::select(
    ID,
    Logistic_AICc, Gompertz_AICc, Baranyi_AICc,
    Logistic_BIC,  Gompertz_BIC,  Baranyi_BIC
  )

three_phase_results <- read_csv("../results/three_phase_model_comparison.csv") %>%
  dplyr::select(
    ID,
    ThreePhase_AICc, ThreePhase_BIC
  )

final_results <- linear_results %>%
  full_join(nonlinear_results, by = "ID") %>%
  full_join(three_phase_results, by = "ID")

write_csv(final_results, "../results/final_model_comparison.csv")

df <- read.csv("../results/final_model_comparison.csv", stringsAsFactors = FALSE)
df[is.na(df)] <- 1e6

model_aicc_cols <- c(
  "Quadratic_AICc", "Cubic_AICc",
  "Logistic_AICc", "Gompertz_AICc", "Baranyi_AICc",
  "ThreePhase_AICc"
)
model_bic_cols <- c(
  "Quadratic_BIC", "Cubic_BIC",
  "Logistic_BIC", "Gompertz_BIC", "Baranyi_BIC",
  "ThreePhase_BIC"
)

best_model_aicc <- apply(df[, model_aicc_cols], 1, function(x) {
  c("Quadratic", "Cubic", "Logistic", "Gompertz", "Baranyi", "ThreePhase")[which.min(x)]
})

best_model_bic <- apply(df[, model_bic_cols], 1, function(x) {
  c("Quadratic", "Cubic", "Logistic", "Gompertz", "Baranyi", "ThreePhase")[which.min(x)]
})

df$Best_Model_AICc <- best_model_aicc
df$Best_Model_BIC  <- best_model_bic

write.csv(df, "../results/best_model_comparison_updated.csv", row.names = FALSE)
message("Final best model selection (AICc & BIC) saved to '../results/best_model_comparison_updated.csv'.")

cat("===== HEAD OF df =====\n")
print(head(df))

best_model_counts_aicc <- table(df$Best_Model_AICc)
message("Best Model Selection Counts (AICc):")
print(best_model_counts_aicc)

best_model_counts_bic <- table(df$Best_Model_BIC)
message("Best Model Selection Counts (BIC):")
print(best_model_counts_bic)

model_onee6_counts_aicc <- sapply(model_aicc_cols, function(col) sum(df[[col]] == 1e6))
message("Each model's AICc is 1e6 for these many rows:")
print(model_onee6_counts_aicc)

model_onee6_counts_bic <- sapply(model_bic_cols, function(col) sum(df[[col]] == 1e6))
message("Each model's BIC is 1e6 for these many rows:")
print(model_onee6_counts_bic)

################################
# Barplot for AICc & BIC Model Selection
################################

model_selection_data <- read_csv("../results/best_model_comparison_updated.csv")
aicc_counts <- table(model_selection_data$Best_Model_AICc)
bic_counts  <- table(model_selection_data$Best_Model_BIC)

model_counts_df <- data.frame(
  Model     = rep(names(aicc_counts), 2),
  Count     = c(as.numeric(aicc_counts), as.numeric(bic_counts)),
  Criterion = rep(c("AICc", "BIC"), each = length(aicc_counts))
)

p <- ggplot(model_counts_df, aes(x = Model, y = Count, fill = Criterion)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Model Selection by AICc and BIC",
       x = "Models", y = "Selection Counts") +
  theme_minimal() +
  scale_fill_manual(values = c("AICc" = "#1f77b4", "BIC" = "#ff7f0e")) +
  geom_text(aes(label = Count), vjust = -0.3, position = position_dodge(0.9))

ggsave("../results/model_selection_comparison.png", plot = p, width = 8, height = 6, dpi = 300)
message("Plot saved successfully to '../results/model_selection_comparison.png'")

############ Additional Example: Plotting Fitted Curves for Each ID ############

growth_data <- read_csv("../data/modified_growth_data.csv")
all_ids <- sort(unique(growth_data$ID))

fit_and_select_best_model_6 <- function(df_sub, n_starts = 50) {
  fit_quad <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2), data = df_sub),
    error = function(e) NULL
  )
  fit_cubic <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df_sub),
    error = function(e) NULL
  )
  fit_logistic <- tryCatch(
    multi_start_fit_one_model(df_sub, "Logistic", n_starts),
    error = function(e) NULL
  )
  fit_gompertz <- tryCatch(
    multi_start_fit_one_model(df_sub, "Gompertz", n_starts),
    error = function(e) NULL
  )
  fit_baranyi <- tryCatch(
    multi_start_fit_one_model(df_sub, "Baranyi", n_starts),
    error = function(e) NULL
  )
  fit_threephase <- tryCatch(
    multi_start_fit_three_phase_model(df_sub, n_starts),
    error = function(e) NULL
  )
  
  if(is.null(fit_quad) || is.null(fit_cubic) || is.null(fit_logistic) ||
     is.null(fit_gompertz) || is.null(fit_baranyi) || is.null(fit_threephase)) {
    return(NULL)
  }
  
  model_list <- list(
    Quadratic = fit_quad,
    Cubic = fit_cubic,
    Logistic = fit_logistic,
    Gompertz = fit_gompertz,
    Baranyi = fit_baranyi,
    ThreePhase = fit_threephase
  )
  aicc_values <- sapply(model_list, AICc)
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  
  t_min <- min(df_sub$Time, na.rm = TRUE)
  t_max <- max(df_sub$Time, na.rm = TRUE)
  new_time <- seq(t_min, t_max, length.out = 200)
  
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
  
  return(list(
    predictions = preds_all,
    best_model = best_model_name
  ))
}

output_dir <- "../results/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for(this_id in all_ids) {
  df_sub <- growth_data %>% filter(ID == this_id)
  
  if(nrow(df_sub) < 8) {
    next
  }
  
  fit_result <- fit_and_select_best_model_6(df_sub, n_starts = 150)
  if(is.null(fit_result)) {
    next
  }
  
  pred_df <- fit_result$predictions
  best_name <- fit_result$best_model
  
  p <- ggplot() + 
    geom_point(data = df_sub, aes(x = Time, y = logPopBio)) + 
    geom_line(data = pred_df, aes(x = Time, y = logPopBio, color = Model), size = 0.8) + 
    labs(
      title = paste0("ID = ", this_id, ", Best Model (AICc) = ", best_name),
      x = "Time",
      y = "log(PopBio)"
    ) + 
    theme_minimal() + 
    scale_color_discrete(name = "Model") + 
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold")
    )
  
  png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
  png(png_filename, width = 2000, height = 1500, res = 200)
  print(p)
  dev.off()
  
  message("Output PNG: ", png_filename, " ; Best model = ", best_name)
}
fit_and_select_best_model_6_partial <- function(df_sub, n_starts = 50) {
  
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
  
  # 将所有模型拟合结果存入列表
  model_list <- list(
    Quadratic  = fit_quad,
    Cubic      = fit_cubic,
    Logistic   = fit_logistic,
    Gompertz   = fit_gompertz,
    Baranyi    = fit_baranyi,
    ThreePhase = fit_threephase
  )
  
  # 只保留拟合成功的模型（非 NULL）
  model_list <- model_list[!sapply(model_list, is.null)]
  
  # 如果没有任何模型拟合成功，则返回 NULL
  if(length(model_list) == 0) {
    return(NULL)
  }
  
  # 计算所有拟合成功模型的 AICc 值，选择 AICc 最小的模型
  aicc_values <- sapply(model_list, AICc)  # 使用 MuMIn::AICc
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  
  # 生成预测网格
  t_min <- min(df_sub$Time, na.rm = TRUE)
  t_max <- max(df_sub$Time, na.rm = TRUE)
  new_time <- seq(t_min, t_max, length.out = 200)
  
  # 对每个拟合成功的模型分别生成预测数据
  preds_list <- lapply(names(model_list), function(mod_name) {
    model_fit <- model_list[[mod_name]]
    data.frame(
      Time = new_time,
      logPopBio = predict(model_fit, newdata = data.frame(Time = new_time)),
      Model = mod_name
    )
  })
  
  preds_all <- do.call(rbind, preds_list)
  
  # 返回包含预测数据和最佳模型名称的列表
  return(list(
    predictions = preds_all,
    best_model  = best_model_name,
    models = model_list
  ))
}

# 设置输出结果的新文件夹路径
output_dir <- "../results/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 对所有 ID 进行循环处理
for(this_id in all_ids) {
  # 取出该 ID 对应的数据子集
  df_sub <- growth_data %>% filter(ID == this_id)
  
  # 跳过数据点少于 10 的子集
  if(nrow(df_sub) < 8) {
    next
  }
  
  # 使用新函数进行 6 模型拟合，允许部分模型失败
  fit_result <- fit_and_select_best_model_6_partial(df_sub, n_starts = 150)
  
  # 如果所有模型都拟合失败，仍然绘制原始数据散点图
  if(is.null(fit_result)) {
    p <- ggplot(df_sub, aes(x = Time, y = logPopBio)) +
      geom_point() +
      labs(
        title = paste0("ID = ", this_id, " (No model fitted successfully)"),
        x = "Time",
        y = "log(PopBio)"
      ) +
      theme_minimal()
    
    png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
    png(png_filename, width = 2000, height = 1500, res = 200)
    print(p)
    dev.off()
    
    message("输出PNG (仅原始数据): ", png_filename)
    next
  }
  
  # 获取预测数据和最佳模型名称
  pred_df   <- fit_result$predictions
  best_name <- fit_result$best_model
  
  # 绘制图像：原始散点图 + 所有拟合成功模型的预测曲线
  p <- ggplot() + 
    geom_point(data = df_sub, aes(x = Time, y = logPopBio)) + 
    geom_line(data = pred_df, aes(x = Time, y = logPopBio, color = Model), size = 0.8) + 
    labs(
      title = paste0("ID = ", this_id, ", Best Model (AICc) = ", best_name),
      x = "Time",
      y = "log(PopBio)"
    ) + 
    theme_minimal() + 
    scale_color_discrete(name = "Model") + 
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold")
    )
  
  # 保存图像为 PNG 文件
  png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
  png(png_filename, width = 2000, height = 1500, res = 200)
  print(p)
  dev.off()  
  
  message("输出PNG: ", png_filename, " ; Best model = ", best_name) 
}

message("全部完成: 现在已有包括三相线性模型在内的 6 模型对比结果及图像保存在 ", output_dir)

message("All tasks completed: Fitted curves and comparisons for 6 models (including three-phase linear model) have been saved.")
