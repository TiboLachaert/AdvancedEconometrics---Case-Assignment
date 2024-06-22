


library(Matrix)

gmm_estimator <- function(y, X, Z) {
  # Ensure the dimensions of y, X, Z are consistent
  T <- dim(y)[1]
  N <- dim(y)[2]
  K <- dim(X)[3]
  R <- dim(Z)[3]
  
  # First difference the variables
  dy <- diff(y, differences = 1)
  dX <- array(0, dim = c(T-1, N, K))
  dZ <- array(0, dim = c(T-1, N, R))
  
  for (k in 1:K) {
    dX[,,k] <- diff(X[,,k], differences = 1)
  }
  
  for (r in 1:R) {
    dZ[,,r] <- diff(Z[,,r], differences = 1)
  }
  
  # Reshape the differenced arrays for matrix operations
  dy_vec <- as.vector(dy)
  dX_mat <- do.call(cbind, lapply(1:N, function(i) dX[,i,]))
  dZ_mat <- do.call(cbind, lapply(1:N, function(i) dZ[,i,]))
  
  # One-step GMM estimation
  W <- diag(ncol(dZ_mat))  # Initial weights matrix (identity)
  beta_one_step <- solve(t(dX_mat) %*% dZ_mat %*% W %*% t(dZ_mat) %*% dX_mat) %*% (t(dX_mat) %*% dZ_mat %*% W %*% t(dZ_mat) %*% dy_vec)
  
  # Residuals from the one-step estimate
  residuals <- dy_vec - dX_mat %*% beta_one_step
  
  # Compute optimal weight matrix
  S <- t(dZ_mat) %*% diag(residuals^2) %*% dZ_mat
  W_opt <- solve(S)
  
  # Two-step GMM estimation
  beta_two_step <- solve(t(dX_mat) %*% dZ_mat %*% W_opt %*% t(dZ_mat) %*% dX_mat) %*% (t(dX_mat) %*% dZ_mat %*% W_opt %*% t(dZ_mat) %*% dy_vec)
  
  # Standard errors, t-stats, and p-values
  V_beta_one_step <- solve(t(dX_mat) %*% dZ_mat %*% W %*% t(dZ_mat) %*% dX_mat)
  se_beta_one_step <- sqrt(diag(V_beta_one_step))
  t_stats_one_step <- beta_one_step / se_beta_one_step
  p_values_one_step <- 2 * pt(-abs(t_stats_one_step), df = T * N - K)
  
  V_beta_two_step <- solve(t(dX_mat) %*% dZ_mat %*% W_opt %*% t(dZ_mat) %*% dX_mat)
  se_beta_two_step <- sqrt(diag(V_beta_two_step))
  t_stats_two_step <- beta_two_step / se_beta_two_step
  p_values_two_step <- 2 * pt(-abs(t_stats_two_step), df = T * N - K)
  
  # Sargan/Hansen test for overidentification
  J_statistic <- t(residuals) %*% dZ_mat %*% W_opt %*% t(dZ_mat) %*% residuals
  p_value_J <- 1 - pchisq(J_statistic, df = R - K)
  
  results <- list(
    one_step = list(
      estimates = beta_one_step,
      standard_errors = se_beta_one_step,
      t_stats = t_stats_one_step,
      p_values = p_values_one_step
    ),
    two_step = list(
      estimates = beta_two_step,
      standard_errors = se_beta_two_step,
      t_stats = t_stats_two_step,
      p_values = p_values_two_step,
      Sargan_Hansen_J = J_statistic,
      p_value_J = p_value_J
    )
  )
  
  return(results)
}



# Example usage:
y <- matrix(rnorm(100), nrow = 10, ncol = 10)
X <- array(rnorm(100 * 3), dim = c(10, 10, 3))
Z <- array(rnorm(100 * 4), dim = c(10, 10, 4))
results <- gmm_estimator(y, X, Z)
print(results)


















################################################
# BALTAGI DATA
################################################

library(tidyverse)
library(readr)
library(fastDummies)
library(matlib)

setwd("D:/users/nicolas.romero/OneDrive - Vlerick Business School/UGent/23-24/econometrics_2/take_home_exam")

data <- read_delim("Data_Baltagi.csv", delim = ";", 
                   escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>% 
  rename("lnC_it" = "ln C_it", 
         "lnP_it" = "ln P_it", 
         "lnPn_it" = "ln Pn_it", 
         "lnY_it" = "ln Y_it") %>%
  group_by(state) %>%
  mutate(lnC_it_lag = lag(lnC_it, n = 1)) %>%
  ungroup() %>%
  na.omit()

# Define y, X for the GMM estimator

## Replicate FE results ----
y_GMM <- data %>% 
  select(state, year, lnC_it) %>%
  pivot_wider(names_from = state, values_from = lnC_it) %>%
  select(-year) %>%
  as.matrix()

X_GMM <- data %>% 
  select(state, year, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  mutate(dummy_cols(data %>% select(year))) %>%
  pivot_longer(cols = -c(state, year)) %>%
  xtabs(data = ., value ~ year + state + name)






