#=====================================#
#     Take home exam 2024 - Part 1    # 
# Tibo Lachaert & Nicolas Romero Diaz #
#=====================================#

# Import functions and packages ----
source("OwnFunctions.R")
library(tidyverse)
library(readr)
library(fastDummies)
library(matlib)
library(openxlsx)

# Case 0: OLS estimation -----

## Import data
data <- read_delim("Data_Baltagi.csv", delim = ";", 
                   escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>% 
  rename("lnC_it" = "ln C_it", 
         "lnP_it" = "ln P_it", 
         "lnPn_it" = "ln Pn_it", 
         "lnY_it" = "ln Y_it") %>%
  group_by(state) %>%
  mutate(intercept = 1) %>%
  mutate(lnC_it_lag = lag(lnC_it, n = 1)) %>%
  mutate(lnC_it_lag2 = lag(lnC_it, n = 2)) %>%
  mutate(lnC_it_lag3 = lag(lnC_it, n = 3)) %>%
  mutate(lnC_it_diff = c(NA, diff(lnC_it, lag = 1))) %>%
  mutate(lnC_it_diff_lag = lag(lnC_it_diff, n = 1)) %>%
  ungroup()

## Replicate OLS results ----
data_OLS <- data %>%
  select(lnC_it, intercept, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  na.omit()

df_OLS_y <- data_OLS %>%
  select(lnC_it) %>%
  as.matrix()

df_OLS_X <- data_OLS %>%  
  select(intercept, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  as.matrix()

OLS <- OLS_own(df_OLS_y, df_OLS_X)
show(OLS[2:5, ])

# Case 1: Fixed Effects (FE) estimation ----

## Replicate FE results ----
data_FE <- data %>%
  select(state, year, lnC_it, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  na.omit()

df_FE_y <- data_FE %>% 
  select(state, year, lnC_it) %>%
  pivot_wider(names_from = state, values_from = lnC_it) %>%
  select(-year) %>%
  as.matrix()

df_FE_X <- data_FE %>% 
  select(state, year, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  mutate(dummy_cols(data_FE %>% select(year))) %>%
  pivot_longer(cols = -c(state, year)) %>%
  xtabs(data = ., value ~ year + state + name)

FE <- FE_own(df_FE_y, df_FE_X)
show(FE[1:4, ]) #Correction term added in function to match t-stats of assignment

FE <- FE_own(df_FE_y, df_FE_X, correction = TRUE)
show(FE[1:4, ]) #Correction term added in function to match t-stats of assignment

## Simulating properties of the FE estimator ----

# Function to generate AR(1) panel data
generate_ar1_panel <- function(rho, N, T) {
  
  # Number of total time steps including burn-in period
  total_T <- T + 25
  
  # Initialize a matrix to hold the panel data
  panel_data <- matrix(0, nrow = N * T, ncol = 3)
  colnames(panel_data) <- c("individual", "time", "value")
  
  for (i in 1:N) {
    # Initialize the AR(1) process for each individual
    x <- numeric(total_T)
    x[1] <- rnorm(1)  # Initial value
    
    # Generate the AR(1) process with burn-in period
    for (t in 2:total_T) {
      x[t] <- rho * x[t - 1] + rnorm(1)
    }
    
    # Remove the burn-in period
    x <- x[-(1:25)]
    
    # Populate the panel data matrix
    panel_data[((i-1) * T + 1):(i * T), ] <- cbind(rep(i, T), 1:T, x)
  }
  
  # Convert the matrix to a dataframe
  panel_df <- as.data.frame(panel_data)
  
  # Convert individual and time columns to factors for panel data representation
  panel_df$individual <- as.factor(panel_df$individual)
  panel_df$time <- as.factor(panel_df$time)
  
  return(panel_df)
}

#Perform Monte Carlo
N_sim <- 100 #Number of simulations
rho_par <- c(0, 0.5, 0.9)
T_par = c(4, 10, 20, 50)
N_par = c(10, 100)

set.seed(11132)
result_FE_MC <- data.frame(rho = double(), T = integer(), N = integer(), 
                        result_rho = double(), result_ese = double(), 
                        result_tse = double(), result_an = double(), 
                        size_an = double(), size_ese = double())
#ese = estimated standard error
#tse = true standard error
#an = analytical standard error

for (rho in rho_par) {
  for (T in T_par) {
    for (N in N_par) {
      print(c(rho, T, N))
      rho_sim  <- c()
      se_sim   <- c()
      
      for (i in 1:N_sim) {
        df_MC   <- generate_ar1_panel(rho, N, T) %>%
          group_by(individual) %>%
          mutate(yit = value) %>%
          mutate(yit_lag1 = lag(yit, n = 1)) %>%
          na.omit()
        
        df_FE_x <- df_MC %>%
          select(individual, time, yit_lag1)  %>%
          pivot_longer(cols = -c(individual, time)) %>%
          mutate(time = as.numeric(time) - 2) %>%
          xtabs(data = ., value ~ time + individual + name)
        
        df_FE_y <- df_MC %>%
          select(individual, time, yit) %>%
          pivot_wider(names_from = individual, values_from = yit) %>%
          select(-time) %>%
          as.matrix()
        
        FE_MC   <- FE_own(df_FE_y, df_FE_x)
        rho_sim <- c(rho_sim, FE_MC[1])
        se_sim  <- c(se_sim, FE_MC[2])
      }
      result_FE_MC <- rbind(result_FE_MC, list(rho = rho, T = T, N = N, 
                                         result_rho = mean(rho_sim), result_ese = mean(se_sim), 
                                         result_tse = sd(rho_sim), result_an = 1/sqrt(T * N * (1 - rho^2)), 
                                         size_an = sum(abs((rho_sim - rho)/(1/sqrt(T * N * (1 - rho^2)))) > 1.96) / N_sim,
                                         size_ese = sum(abs((rho_sim - rho)/se_sim) > 1.96) / N_sim))
    }
  }
}

result_FE_MC$asymptotic_assumption <- result_FE_MC$N / result_FE_MC$T^3 
result_FE_MC$nickel_bias           <- - (1+result_FE_MC$rho)/result_FE_MC$T
result_FE_MC$nickel_bias_corrected <- result_FE_MC$rho + result_FE_MC$nickel_bias

## Simulation methods for bias correction and inference ----
#See jackknife_bias_correction.R

### Half-panel jackknife bias-corrected FE estimator ----

### Bootstrap ----

# Case 2: Generalized Method of Moments (GMM) estimation ----

## Test for one setting ----
rho <- 0.5
N <- 10
T <- 20

set.seed(11132)
df_GMM <- generate_ar1_panel(rho, N, T) %>% 
  group_by(individual) %>%
  mutate(yit = value) %>%
  mutate(yit_lag1 = lag(yit, n = 1)) %>%
  mutate(yit_lag2 = lag(yit, n = 2)) %>%
  mutate(yit_lag3 = lag(yit, n = 3)) %>%
  mutate(yit_diff = yit - lag(yit)) %>%
  mutate(yit_diff_lag = lag(yit_diff)) %>%
  ungroup()

df_GMM_IV1 <- df_GMM %>%
  select(individual, time, yit_diff, yit_diff_lag, yit_lag2) %>%
  na.omit()

df_GMM_IV1_y <- df_GMM_IV1 %>% 
  select(individual, time, yit_diff) %>%
  pivot_wider(names_from = individual, values_from = yit_diff) %>%
  select(-time) %>%
  as.matrix()

df_GMM_IV1_X <- df_GMM_IV1 %>% 
  select(individual, time, yit_diff_lag) %>%
  pivot_longer(cols = -c(individual, time)) %>%
  mutate(time = as.numeric(time) - 2) %>%
  xtabs(data = ., value ~ time + individual + name)

df_GMM_IV1_Z <- df_GMM_IV1 %>% 
  select(individual, time, yit_lag2) %>%
  pivot_longer(cols = -c(individual, time)) %>%
  mutate(time = as.numeric(time) - 2) %>%
  xtabs(data = ., value ~ time + individual + name)

df_GMM_IV2 <- df_GMM %>%
  select(individual, time, yit_diff, yit_diff_lag, yit_lag2, yit_lag3) %>%
  na.omit()

df_GMM_IV2_y <- df_GMM_IV2 %>% 
  select(individual, time, yit_diff) %>%
  pivot_wider(names_from = individual, values_from = yit_diff) %>%
  select(-time) %>%
  as.matrix()

df_GMM_IV2_X <- df_GMM_IV2 %>% 
  select(individual, time, yit_diff_lag) %>%
  pivot_longer(cols = -c(individual, time)) %>%
  mutate(time = as.numeric(time) - 3) %>%
  xtabs(data = ., value ~ time + individual + name)

df_GMM_IV2_Z <- df_GMM_IV2 %>% 
  select(individual, time, yit_lag2, yit_lag3) %>%
  pivot_longer(cols = -c(individual, time)) %>%
  mutate(time = as.numeric(time) - 3) %>%
  xtabs(data = ., value ~ time + individual + name)

GMM <- GMM_own(df_GMM_IV1_y, df_GMM_IV1_X, df_GMM_IV1_Z)
show(GMM)

GMM <- GMM_own(df_GMM_IV2_y, df_GMM_IV2_X, df_GMM_IV2_Z)
show(GMM)

## Monte Carlo ----

### With one instrument ----

onestep_MC <- function(N_sim, rho, N, T){
  results <- matrix(data = 0, nrow = N_sim, ncol = 5)
  colnames(results) <- c("coefs", "stdvs", "tstats", "pvals", "Jstat")
  
  for (i in 1:N_sim){
    df_GMM <- generate_ar1_panel(rho, N, T) %>% 
      group_by(individual) %>%
      mutate(yit = value) %>%
      mutate(yit_lag1 = lag(yit, n = 1)) %>%
      mutate(yit_lag2 = lag(yit, n = 2)) %>%
      mutate(yit_lag3 = lag(yit, n = 3)) %>%
      mutate(yit_diff = yit - lag(yit)) %>%
      mutate(yit_diff_lag = lag(yit_diff)) %>%
      ungroup()
    
    df_GMM_IV1 <- df_GMM %>%
      select(individual, time, yit_diff, yit_diff_lag, yit_lag2) %>%
      na.omit()
    
    df_GMM_IV1_y <- df_GMM_IV1 %>% 
      select(individual, time, yit_diff) %>%
      pivot_wider(names_from = individual, values_from = yit_diff) %>%
      select(-time) %>%
      as.matrix()
    
    df_GMM_IV1_X <- df_GMM_IV1 %>% 
      select(individual, time, yit_diff_lag) %>%
      pivot_longer(cols = -c(individual, time)) %>%
      mutate(time = as.numeric(time) - 2) %>%
      xtabs(data = ., value ~ time + individual + name)
    
    df_GMM_IV1_Z <- df_GMM_IV1 %>% 
      select(individual, time, yit_lag2) %>%
      pivot_longer(cols = -c(individual, time)) %>%
      mutate(time = as.numeric(time) - 2) %>%
      xtabs(data = ., value ~ time + individual + name)
    
    GMM <- GMM_own(df_GMM_IV1_y, df_GMM_IV1_X, df_GMM_IV1_Z)
    
    results[i, ] <- GMM$one_step
  }
  return(c(mean(results[, 1]), mean(results[, 2]), sd(results[, 1]), sum(abs((results[, 1] - rho)/results[, 2]) > 1.96)/N_sim))
}

set.seed(11132)
result_onestep_MC <- data.frame(rho = double(), T = integer(), N = integer(), 
                        result_rho = double(), result_ese = double(), 
                        result_tse = double(), size = double())

for (rho in rho_par) {
  for (T in T_par) {
    for (N in N_par) {
      print(c(rho, T, N))
      result <- onestep_MC(N_sim, rho, N, T)
      
      result_onestep_MC <- rbind(result_onestep_MC, list(rho = rho, T = T, N = N, 
                                         result = result[1], result_ese = result[2], 
                                         result_tse = result[3], size = result[4]))
    }
  }
}

### With two instruments ----

twostep_MC <- function(N_sim, rho, N, T){
  results <- matrix(data = 0, nrow = N_sim, ncol = 5)
  colnames(results) <- c("coefs", "stdvs", "tstats", "pvals", "Jstat")
  
  for (i in 1:N_sim){
    df_GMM <- generate_ar1_panel(rho, N, T) %>% 
      group_by(individual) %>%
      mutate(yit = value) %>%
      mutate(yit_lag1 = lag(yit, n = 1)) %>%
      mutate(yit_lag2 = lag(yit, n = 2)) %>%
      mutate(yit_lag3 = lag(yit, n = 3)) %>%
      mutate(yit_diff = yit - lag(yit)) %>%
      mutate(yit_diff_lag = lag(yit_diff)) %>%
      ungroup()
    
    df_GMM_IV2 <- df_GMM %>%
      select(individual, time, yit_diff, yit_diff_lag, yit_lag2, yit_lag3) %>%
      na.omit()
    
    df_GMM_IV2_y <- df_GMM_IV2 %>% 
      select(individual, time, yit_diff) %>%
      pivot_wider(names_from = individual, values_from = yit_diff) %>%
      select(-time) %>%
      as.matrix()
    
    df_GMM_IV2_X <- df_GMM_IV2 %>% 
      select(individual, time, yit_diff_lag) %>%
      pivot_longer(cols = -c(individual, time)) %>%
      mutate(time = as.numeric(time) - 3) %>%
      xtabs(data = ., value ~ time + individual + name)
    
    df_GMM_IV2_Z <- df_GMM_IV2 %>% 
      select(individual, time, yit_lag2, yit_lag3) %>%
      pivot_longer(cols = -c(individual, time)) %>%
      mutate(time = as.numeric(time) - 3) %>%
      xtabs(data = ., value ~ time + individual + name)
    
    GMM <- GMM_own(df_GMM_IV2_y, df_GMM_IV2_X, df_GMM_IV2_Z)
    
    results[i, ] <- GMM$two_step
  }
  return(c(mean(results[, 1]), mean(results[, 2]), sd(results[, 1]), sum(abs((results[, 1] - rho)/results[, 2]) > 1.96)/N_sim))
}

set.seed(11132)
result_twostep_MC <- data.frame(rho = double(), T = integer(), N = integer(), 
                                result_rho = double(), result_ese = double(), 
                                result_tse = double(), size = double())

for (rho in rho_par) {
  for (T in T_par) {
    for (N in N_par) {
      print(c(rho, T, N))
      result <- twostep_MC(N_sim, rho, N, T)
      
      result_twostep_MC <- rbind(result_twostep_MC, list(rho = rho, T = T, N = N, 
                                                         result = result[1], result_ese = result[2], 
                                                         result_tse = result[3], size = result[4]))
    }
  }
}

# Export results to Excel ----
names <- list('result_FE_MC' = result_FE_MC, 'result_onestep_MC' = result_onestep_MC, 'result_twostep_MC' = result_twostep_MC)
write.xlsx(names, file = 'Case assignment - Results.xlsx')