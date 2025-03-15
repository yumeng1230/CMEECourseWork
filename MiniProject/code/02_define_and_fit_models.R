#############################################
# 02_define_and_fit_models.R
# 1) Read cleaned data
# 2) Define nonlinear models (Logistic, Gompertz, Baranyi, Three-phase)
# 3) Define multi-start nlsLM fitting functions
# 4) Fit these nonlinear models in a standard way (if desired), or
#    Provide them for later usage by other scripts
#############################################


library(dplyr)
library(readr)
library(minpack.lm)  # for nlsLM
library(MuMIn)       # for AICc
library(broom)       # for glance, tidy

# 1) Read the cleaned dataset
cleaned_data <- read_csv("../data/modified_growth_data.csv")

# 2) Define nonlinear model formulas
logistic_model <- function(t, r_max, K, N_0) {
  log(N_0 * K * exp(r_max * t) / (K + N_0 * (exp(r_max * t) - 1)))
}

gompertz_model <- function(t, r_max, K, N_0, t_lag) {
  N_0 + (K - N_0) * exp(
    -exp(r_max * exp(1) * (t_lag - t) / ((K - N_0) * log(10)) + 1)
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

three_phase_linear_model <- function(t, N0, Nmax, tLAG, tMAX) {
  ifelse(t <= tLAG,
         N0,
         ifelse(t >= tMAX,
                Nmax,
                N0 + (Nmax - N0) / (tMAX - tLAG) * (t - tLAG)
         )
  )
}

# 3) Multi-start fitting approach
compute_max_slope <- function(df) {
  df2 <- df %>% arrange(Time)
  slopes <- diff(df2$logPopBio) / diff(df2$Time)
  slopes <- slopes[!is.na(slopes) & !is.infinite(slopes)]
  if (length(slopes) == 0) return(0.1)
  max(slopes)
}

set.seed(123)
generate_random_inits <- function(n_starts, df) {
  smax <- compute_max_slope(df)
  if (smax < 0.01) smax <- 0.01
  rmax_min <- 0.01
  rmax_max <- min(smax, 5)
  K_min <- min(df$logPopBio, na.rm=TRUE)
  K_max <- max(df$logPopBio, na.rm=TRUE)
  N0_min <- K_min
  N0_max <- K_max
  tlag_min <- 0
  tlag_max <- max(df$Time, na.rm=TRUE)
  tibble(
    r_max = runif(n_starts, rmax_min*1.01, rmax_max*0.99),
    K     = runif(n_starts, K_min*1.01,    K_max*0.99),
    N_0   = runif(n_starts, N0_min*1.01,   N0_max*0.99),
    t_lag = runif(n_starts, tlag_min*1.01, tlag_max*0.99)
  )
}

# For Logistic/Gompertz/Baranyi
multi_start_fit_one_model <- function(df, model_name, n_starts=150) {
  if (nrow(df) < 8) return(NULL)
  form_list <- list(
    Logistic = logPopBio ~ logistic_model(Time, r_max, K, N_0),
    Gompertz = logPopBio ~ gompertz_model(Time, r_max, K, N_0, t_lag),
    Baranyi  = logPopBio ~ baranyi_model(Time, r_max, K, N_0, t_lag)
  )
  formula_used <- form_list[[model_name]]
  param_names <- if (model_name=="Logistic") {
    c("r_max","K","N_0")
  } else {
    c("r_max","K","N_0","t_lag")
  }
  inits <- generate_random_inits(n_starts, df)
  fits_list <- list()
  for(i in seq_len(n_starts)) {
    start_i <- as.list(inits[i, param_names, drop=FALSE])
    res <- tryCatch(
      nlsLM(formula_used, data=df, start=start_i),
      error=function(e) NULL
    )
    if(!is.null(res)) {
      fits_list[[length(fits_list)+1]] <- res
    }
  }
  if(length(fits_list)==0) return(NULL)
  # pick best by AICc
  best_fit <- NULL
  best_aicc <- Inf
  for(f in fits_list) {
    cur_aicc <- AICc(f)
    if(cur_aicc < best_aicc) {
      best_aicc <- cur_aicc
      best_fit <- f
    }
  }
  best_fit
}

# For Three-phase linear model
multi_start_fit_three_phase_model <- function(df, n_starts=150) {
  if(nrow(df)<8) return(NULL)
  form_3p <- logPopBio ~ three_phase_linear_model(Time, N0, Nmax, tLAG, tMAX)
  t_min <- min(df$Time)
  t_max <- max(df$Time)
  log_min <- min(df$logPopBio)
  log_max <- max(df$logPopBio)
  init_mat <- replicate(n_starts, {
    tt <- sort(runif(2, t_min, t_max))
    nn <- sort(runif(2, log_min, log_max))
    c(N0=nn[1], Nmax=nn[2], tLAG=tt[1], tMAX=tt[2])
  })
  fits_list <- list()
  for(i in seq_len(n_starts)) {
    start_i <- as.list(init_mat[,i])
    res <- tryCatch(
      nlsLM(form_3p, data=df, start=start_i),
      error=function(e) NULL
    )
    if(!is.null(res)) {
      fits_list[[length(fits_list)+1]] <- res
    }
  }
  if(length(fits_list)==0) return(NULL)
  # pick best by AICc
  best_fit <- NULL
  best_aicc <- Inf
  for(f in fits_list) {
    cur_aicc <- AICc(f)
    if(cur_aicc < best_aicc) {
      best_aicc <- cur_aicc
      best_fit <- f
    }
  }
  best_fit
}

message("All nonlinear models and multi-start fitting functions are defined.")

