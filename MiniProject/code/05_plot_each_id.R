#############################################

# 05_plot_each_id_partial.R

# 1) Fit 6 models (Quadratic, Cubic, Logistic, Gompertz, Baranyi, Three-phase)

#    allowing partial success

# 2) If all fail, plot only raw points

# 3) Otherwise, plot data points + all successful fits

#############################################



library(dplyr)

library(ggplot2)

library(readr)

# 1) Read data

growth_data <- read_csv("../data/modified_growth_data.csv")

all_ids <- unique(growth_data$ID)

# 2) Source the existing model-fitting functions

source("02_define_and_fit_models.R")

# A function that attempts 6 models; only keeps successful fits

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
  
  # Put all fits in a list
  
  model_list <- list(
    
    Quadratic  = fit_quad,
    
    Cubic      = fit_cubic,
    
    Logistic   = fit_logistic,
    
    Gompertz   = fit_gompertz,
    
    Baranyi    = fit_baranyi,
    
    ThreePhase = fit_threephase
    
  )
  
  # Keep only successful (non-NULL)
  
  model_list <- model_list[!sapply(model_list, is.null)]
  
  # If none succeeded, return NULL
  
  if(length(model_list) == 0) {
    
    return(NULL)
    
  }
  
  # Compute AICc for each successful model & pick best
  
  # (Quadratic/Cubic are lm, others are nls -> all support AICc)
  
  aicc_values <- sapply(model_list, AICc) 
  
  best_model_name <- names(aicc_values)[which.min(aicc_values)]
  
  # Make a time grid
  
  t_min <- min(df_sub$Time, na.rm=TRUE)
  
  t_max <- max(df_sub$Time, na.rm=TRUE)
  
  new_time <- seq(t_min, t_max, length.out = 200)
  
  # Predict from each successful model
  
  preds_list <- lapply(names(model_list), function(mod_name) {
    
    fit_obj <- model_list[[mod_name]]
    
    data.frame(
      
      Time = new_time,
      
      logPopBio = predict(fit_obj, newdata=data.frame(Time=new_time)),
      
      Model = mod_name
      
    )
    
  })
  
  preds_all <- do.call(rbind, preds_list)
  
  # Return a list of predictions + best model name
  
  list(
    
    predictions = preds_all,
    
    best_model  = best_model_name,
    
    models      = model_list
    
  )
  
}

# 3) Create an output directory for the ID plots

out_dir <- "../results/"

if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

for(this_id in all_ids) {
  
  df_sub <- growth_data %>% filter(ID == this_id)
  
  # Skip if too few data points
  
  if(nrow(df_sub) < 8){
    
    next
    
  }
  
  # Attempt partial-fitting approach
  
  fit_result <- fit_and_select_best_model_6_partial(df_sub, n_starts=150)
  
  # If no model succeeded, plot only raw points
  
  if(is.null(fit_result)) {
    
    p <- ggplot(df_sub, aes(x=Time, y=logPopBio)) +
      
      geom_point() +
      
      labs(
        
        title = paste0("ID = ", this_id, " (No model fitted successfully)"),
        
        x = "Time",
        
        y = "log(PopBio)"
        
      ) +
      
      theme_minimal()
    
    png_file <- paste0(out_dir, "6models_ID_", this_id, ".png")
    
    png(png_file, width=2000, height=1500, res=200)
    
    print(p)
    
    dev.off()
    
    message("Output PNG (only data): ", png_file)
    
    next
    
  }
  
  # Otherwise, we have at least one successful model
  
  pred_df <- fit_result$predictions
  
  best_name <- fit_result$best_model
  
  # Plot data + lines for all successful models
  
  p <- ggplot() +
    
    geom_point(data=df_sub, aes(x=Time, y=logPopBio)) +
    
    geom_line(data=pred_df, aes(x=Time, y=logPopBio, color=Model), size=0.8) +
    
    labs(
      
      title = paste0("ID = ", this_id, ", Best Model (AICc) = ", best_name),
      
      x="Time", y="log(PopBio)"
      
    ) +
    
    theme_minimal() +
    
    scale_color_discrete(name="Model") +
    
    theme(
      
      legend.text=element_text(size=14),
      
      legend.title=element_text(size=16, face="bold")
      
    )
  
  # Save the figure
  
  png_file <- paste0(out_dir, "6models_ID_", this_id, ".png")
  
  png(png_file, width=2000, height=1500, res=200)
  
  print(p)
  
  dev.off()
  
  message("Output PNG: ", png_file, " ; Best model = ", best_name)
  
}

message("All ID plots completed. Check ../results/ for output.")

