# New version: Even if some models fail to converge, still plot the results
# and save them into a new folder

source("DataClean.R")
source("Main.R")
source("Parameter.R")
source("Parameter_visualization.R")
fit_and_select_best_model_6_partial <- function(df_sub, n_starts = 50) {
  # 1) Quadratic model
  fit_quad <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2), data = df_sub),
    error = function(e) NULL
  )
  # 2) Cubic model
  fit_cubic <- tryCatch(
    lm(logPopBio ~ Time + I(Time^2) + I(Time^3), data = df_sub),
    error = function(e) NULL
  )
  # 3) Logistic model
  fit_logistic <- tryCatch(
    multi_start_fit_one_model(df_sub, "Logistic", n_starts),
    error = function(e) NULL
  )
  # 4) Modified Gompertz model
  fit_gompertz <- tryCatch(
    multi_start_fit_one_model(df_sub, "Gompertz", n_starts),
    error = function(e) NULL
  )
  # 5) Baranyi model
  fit_baranyi <- tryCatch(
    multi_start_fit_one_model(df_sub, "Baranyi", n_starts),
    error = function(e) NULL
  )
  # 6) Three-phase linear model
  fit_threephase <- tryCatch(
    multi_start_fit_three_phase_model(df_sub, n_starts),
    error = function(e) NULL
  )
  # Combine all fitted models into a list
  model_list <- list(
    Quadratic  = fit_quad,
    Cubic      = fit_cubic,
    Logistic   = fit_logistic,
    Gompertz   = fit_gompertz,
    Baranyi    = fit_baranyi,
    ThreePhase = fit_threephase
  )
  # Keep only those models which successfully converged (not NULL)
  model_list <- model_list[!sapply(model_list, is.null)]
  # If no model has converged, return NULL
  if (length(model_list) == 0) {
    return(NULL)
  }
  # Calculate AICc values for the successfully fitted models and select the best model
  aicc_values <- sapply(model_list, AICc)  # from MuMIn::AICc
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  # Generate a time sequence for predictions
  t_min <- min(df_sub$Time, na.rm = TRUE)
  t_max <- max(df_sub$Time, na.rm = TRUE)
  new_time <- seq(t_min, t_max, length.out = 200)
  # Create prediction data for each successfully fitted model
  preds_list <- lapply(names(model_list), function(mod_name) {
    model_fit <- model_list[[mod_name]]
    data.frame(
      Time = new_time,
      logPopBio = predict(model_fit, newdata = data.frame(Time = new_time)),
      Model = mod_name
    )
  })
  preds_all <- do.call(rbind, preds_list)
  # Return the prediction data and the name of the best model
  return(list(
    predictions = preds_all,
    best_model  = best_model_name,
    models = model_list
  ))
}

# Specify the output directory for results
output_dir <- "../results_partial/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Loop through all IDs
for(this_id in all_ids) {
  # Subset data for the current ID
  df_sub <- growth_data %>% filter(ID == this_id)
  # Skip subsets with fewer than 8 data points
  if(nrow(df_sub) < 8) {
    next
  }
  # Use the new function to fit 6 models; partial failures are allowed
  fit_result <- fit_and_select_best_model_6_partial(df_sub, n_starts = 150)
  # If no models converged, still plot the raw data
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
    message("Output PNG (only raw data): ", png_filename)
    next
  }
  # Extract the prediction dataframe and best-model name
  pred_df   <- fit_result$predictions
  best_name <- fit_result$best_model
  # Create a plot: raw data points + predictions from all converged models
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
  # Save the plot as a PNG
  png_filename <- paste0(output_dir, "6models_ID_", this_id, ".png")
  png(png_filename, width = 2000, height = 1500, res = 200)
  print(p)
  dev.off()  
  message("Output PNG: ", png_filename, " ; Best model = ", best_name) 
}

message("All done: The 6-model comparison (including the three-phase linear) results and plots are saved in ", output_dir)