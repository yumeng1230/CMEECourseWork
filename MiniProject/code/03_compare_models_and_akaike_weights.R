#############################################

# 03_compare_models_and_akaike_weights.R

# 1) Compare linear vs. nonlinear (Logistic, Gompertz, Baranyi, Three-phase)

# 2) Select best models by AICc/BIC

# 3) Compute Akaike weights

# 4) Produce comparison plots

#############################################


library(dplyr)

library(readr)

library(MuMIn)

library(ggplot2)

library(broom)

# 1) Read data and results from 01 & 02

cleaned_data <- read_csv("../data/modified_growth_data.csv")

res_linear   <- read_csv("../results/model_comparison_linear.csv")

# Source the model definitions (02) if you need re-fitting:

source("02_define_and_fit_models.R")

# 2) Strictly fit the 4 nonlinear models for each ID

fit_nonlinear_strict <- function(df, n_starts=150) {
  
  if(nrow(df) < 8) return(NULL)
  
  m_log   <- multi_start_fit_one_model(df,"Logistic",  n_starts)
  
  m_gomp  <- multi_start_fit_one_model(df,"Gompertz",  n_starts)
  
  m_baran <- multi_start_fit_one_model(df,"Baranyi",   n_starts)
  
  m_3p    <- multi_start_fit_three_phase_model(df,     n_starts)
  
  info <- function(f) {
    
    if(is.null(f)) return(c(NA,NA,NA))
    
    c(AICc(f), BIC(f), glance(f)$r.squared)
    
  }
  
  tibble(
    
    ID = unique(df$ID),
    
    Logistic_AICc = info(m_log)[1],
    
    Logistic_BIC  = info(m_log)[2],
    
    Logistic_R2   = info(m_log)[3],
    
    Gompertz_AICc = info(m_gomp)[1],
    
    Gompertz_BIC  = info(m_gomp)[2],
    
    Gompertz_R2   = info(m_gomp)[3],
    
    Baranyi_AICc  = info(m_baran)[1],
    
    Baranyi_BIC   = info(m_baran)[2],
    
    Baranyi_R2    = info(m_baran)[3],
    
    ThreePhase_AICc = info(m_3p)[1],
    
    ThreePhase_BIC  = info(m_3p)[2],
    
    ThreePhase_R2   = info(m_3p)[3]
    
  )
  
}

results_nonlinear <- cleaned_data %>%
  
  group_by(ID) %>%
  
  group_split() %>%
  
  lapply(fit_nonlinear_strict, n_starts=150) %>%
  
  bind_rows()

write_csv(results_nonlinear, "../results/model_comparison_nonlinear.csv")

message("Nonlinear model comparison saved to ../results/model_comparison_nonlinear.csv")

# 3) Merge linear & nonlinear, select best model by AICc/BIC

final_res <- full_join(res_linear, results_nonlinear, by="ID")

final_res[is.na(final_res)] <- 1e6  # large penalty for failed fits

model_aicc_cols <- c("Quadratic_AICc","Cubic_AICc",
                     
                     "Logistic_AICc","Gompertz_AICc",
                     
                     "Baranyi_AICc","ThreePhase_AICc")

model_bic_cols  <- c("Quadratic_BIC","Cubic_BIC",
                     
                     "Logistic_BIC","Gompertz_BIC",
                     
                     "Baranyi_BIC","ThreePhase_BIC")

model_names <- c("Quadratic","Cubic","Logistic","Gompertz","Baranyi","ThreePhase")

best_model_aicc <- apply(final_res[, model_aicc_cols], 1, function(x){
  
  model_names[which.min(x)]
  
})

best_model_bic <- apply(final_res[, model_bic_cols], 1, function(x){
  
  model_names[which.min(x)]
  
})

final_res$Best_Model_AICc <- best_model_aicc

final_res$Best_Model_BIC  <- best_model_bic

write_csv(final_res, "../results/best_model_comparison_updated.csv")

message("Combined results written to ../results/best_model_comparison_updated.csv")

# 4) Plot counts of best models

aicc_counts <- table(final_res$Best_Model_AICc)

bic_counts  <- table(final_res$Best_Model_BIC)

df_counts <- data.frame(
  
  Model=rep(names(aicc_counts),2),
  
  Count=c(as.numeric(aicc_counts), as.numeric(bic_counts)),
  
  Criterion=rep(c("AICc","BIC"), each=length(aicc_counts))
  
)

p1 <- ggplot(df_counts, aes(x=Model, y=Count, fill=Criterion)) +
  
  geom_bar(stat="identity", position="dodge") +
  
  geom_text(aes(label=Count), position=position_dodge(0.9), vjust=-0.3) +
  
  scale_fill_manual(values=c("#1f77b4","#ff7f0e")) +
  
  theme_minimal() +
  
  labs(
    
    title="Best Model Selection Counts (AICc vs BIC)",
    
    x="Model", y="Number of IDs"
    
  )

ggsave("../results/model_selection_comparison.png", p1, width=8, height=6, dpi=300)

message("Bar plot saved: ../results/model_selection_comparison.png")

# 5) Akaike weights

calc_akaike_weights <- function(aic_vec){
  
  delta <- aic_vec - min(aic_vec)
  
  rel_lik <- exp(-0.5 * delta)
  
  rel_lik / sum(rel_lik)
  
}

final_res2 <- final_res

final_res2[final_res2==1e6] <- Inf

weights_matrix <- t(apply(final_res2[, model_aicc_cols], 1, calc_akaike_weights))

colnames(weights_matrix) <- paste0(model_names, "_w")

df_with_weights <- cbind(final_res2, weights_matrix)

write_csv(df_with_weights, "../results/final_model_comparison_with_akaike_weights.csv")

# The best model by weight:

df_with_weights$Best_Model_Weights <- apply(
  
  df_with_weights[, paste0(model_names,"_w")],
  
  1, function(x) model_names[which.max(x)]
  
)

# Plot Akaike weight distribution

library(tidyr)

df_w_long <- df_with_weights %>%
  dplyr::select(ID, all_of(paste0(model_names,"_w"))) %>%
  pivot_longer(-ID, names_to="Model", values_to="Weight")


p2 <- ggplot(df_w_long, aes(x=Model, y=Weight, fill=Model)) +
  
  geom_boxplot(outlier.shape=NA) +
  
  geom_jitter(width=0.2, alpha=0.5, size=1.8) +
  
  theme_minimal() +
  
  theme(legend.position="right")+
  labs(
    title = "Akaike Weight Distribution (6 Models)",
    x = "Model", 
    y = "Akaike Weight",
    fill = "Model"
  )

ggsave("../results/akaike_weights_boxplot.png", p2, width=8, height=6, dpi=300)

message("Akaike weights boxplot saved: ../results/akaike_weights_boxplot.png")

