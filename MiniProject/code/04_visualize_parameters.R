#############################################

# 04_visualize_parameters.R

# 1) Extract parameters from each nonlinear model

# 2) Summarize distribution via boxplots

#############################################

library(dplyr)

library(readr)

library(minpack.lm)

library(MuMIn)

library(broom)

library(ggplot2)

source("02_define_and_fit_models.R")

# 1) Read cleaned data

data_for_params <- read_csv("../data/modified_growth_data.csv")

# Define a function to strictly fit each of the 4 nonlinear models

extract_params_strict <- function(df, n_starts=150) {
  
  if(nrow(df) < 8) return(NULL)
  
  f_log <- multi_start_fit_one_model(df, "Logistic", n_starts)
  
  f_gom <- multi_start_fit_one_model(df, "Gompertz", n_starts)
  
  f_bar <- multi_start_fit_one_model(df, "Baranyi",  n_starts)
  
  f_3ph <- multi_start_fit_three_phase_model(df,     n_starts)
  
  get_tidy <- function(m, mod_name) {
    
    if(is.null(m)) return(NULL)
    
    out <- tidy(m)
    
    out$Model <- mod_name
    
    out
    
  }
  
  out_df <- bind_rows(
    
    get_tidy(f_log, "Logistic"),
    
    get_tidy(f_gom, "Gompertz"),
    
    get_tidy(f_bar, "Baranyi"),
    
    get_tidy(f_3ph, "ThreePhase")
    
  )
  
  if(nrow(out_df)==0) return(NULL)
  
  out_df$ID <- unique(df$ID)
  
  out_df
  
}

# 2) Extract parameters for each ID

params_all <- data_for_params %>%
  
  group_by(ID) %>%
  
  group_split() %>%
  
  lapply(extract_params_strict, n_starts=150) %>%
  
  bind_rows()

write_csv(params_all, "../results/nonlinear_model_parameters.csv")

message("Nonlinear model parameters extracted: ../results/nonlinear_model_parameters.csv")

# 3) Create boxplots

param_df <- read_csv("../results/nonlinear_model_parameters.csv") %>%
  
  mutate(
    
    term = case_when(
      
      term %in% c("N0","N_0")     ~ "N0",
      
      term %in% c("tLAG","t_lag") ~ "t_lag",
      
      TRUE ~ term
      
    )
    
  )

unique_terms <- sort(unique(param_df$term))

my_theme <- theme_bw() +
  
  theme(
    
    text=element_text(size=14),
    
    plot.title=element_text(face="bold", size=16, hjust=0.5),
    
    axis.title=element_text(face="bold", size=14),
    
    axis.text=element_text(size=12),
    
    legend.position="none"
    
  )

for(p_name in unique_terms) {
  
  df_sub <- param_df %>% filter(term==p_name)
  
  p <- ggplot(df_sub, aes(x=Model, y=estimate, fill=Model)) +
    
    geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.6) +
    
    geom_jitter(width=0.2, size=2, alpha=0.6) +
    
    scale_fill_brewer(palette="Set2") +
    
    labs(
      
      title=paste("Parameter", p_name, "Comparison Across Models"),
      
      x="Model", y="Estimated Value"
      
    ) +
    
    my_theme
  
  # Manually adjust y-limits for certain parameters:
  
  if(p_name=="t_lag"){
    
    p <- p + coord_cartesian(ylim=c(-10,100))
    
  } else if(p_name=="r_max"){
    
    p <- p + coord_cartesian(ylim=c(-10,10))
    
  } else if(p_name=="N0"){
    
    p <- p + coord_cartesian(ylim=c(-50,100))
    
  } else if(p_name=="K"){
    
    p <- p + coord_cartesian(ylim=c(-50,100))
    
  }
  
  fig_path <- paste0("../results/nonlinear_parameter_comparison_", p_name, ".png")
  
  ggsave(fig_path, p, width=7, height=5, dpi=300)
  
  message("Parameter boxplot for '", p_name, "' saved: ", fig_path)
  
}

