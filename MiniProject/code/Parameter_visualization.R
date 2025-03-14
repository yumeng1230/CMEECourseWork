##############################
# visualize_nonlinear_parameters_improved_zoom.R
# Compare parameter estimates for nonlinear models individually 
# (merging parameters, beautifying plots + zooming for t_lag)
##############################

rm(list = ls())
library(dplyr)
library(ggplot2)
library(readr)

# 1. Read the parameter data
params <- read_csv("../results/model_parameters_extracted.csv")

# 2. Keep only the nonlinear models
nonlinear_params <- params %>%
  filter(Model %in% c("Logistic", "Gompertz", "Baranyi", "ThreePhase"))

# 3. Standardize parameter names
nonlinear_params <- nonlinear_params %>%
  mutate(term = case_when(
    term %in% c("N0", "N_0") ~ "N0",
    term %in% c("tLAG", "t_lag") ~ "t_lag",
    TRUE ~ term
  ))

# 4. List all parameter names
unique_terms <- sort(unique(nonlinear_params$term))

# 5. Define a custom theme for the plots
custom_theme <- theme_bw() +
  theme(
    text = element_text(family = "sans", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# 6. Generate comparison plots for each parameter
output_dir <- "../results/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (param in unique_terms) {
  df_param <- nonlinear_params %>% filter(term == param)
  
  p <- ggplot(df_param, aes(x = Model, y = estimate, fill = Model)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste("Comparison of parameter", param, "estimates in nonlinear models"),
      x = "Model",
      y = "Estimate"
    ) +
    custom_theme
  
  # If the parameter is t_lag, zoom the y-axis to focus around 0
  # For example, here the limits are set to [-10, 100]; adjust as needed
  if (param == "t_lag") {
    p <- p + coord_cartesian(ylim = c(-10, 100))
  }
  if (param == "r_max") {
    p <- p + coord_cartesian(ylim = c(-10, 10))
  }
  if (param == "N0") {
    p <- p + coord_cartesian(ylim = c(-50, 100))
  }
  if (param == "K") {
    p <- p + coord_cartesian(ylim = c(-50, 100))
  }
  
  # Save the plot
  filename <- paste0(output_dir, "nonlinear_parameter_comparison_", param, ".png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  message("Comparison plot for parameter ", param, " has been saved to: ", filename)
}
