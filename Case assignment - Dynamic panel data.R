# Part 0: OLS estimation -----

## Import data and packages ----
source("OwnFunctions.R")
library(tidyverse)
library(readr)
library(fastDummies)
library(matlib)

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

## Replicate OLS results ----
df_OLS_y <- data %>%
  select(lnC_it) %>%
  as.matrix()

df_OLS_X <- data %>%  
  mutate(intercept = 1) %>%
  select(intercept, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  as.matrix()

OLS <- OLS_own(df_OLS_y, df_OLS_X)
show(OLS[2:5, ])

# Part 1: Fixed Effects (FE) estimation ----

## Replicate FE results ----
df_FE_y <- data %>% 
  select(state, year, lnC_it) %>%
  pivot_wider(names_from = state, values_from = lnC_it) %>%
  select(-year) %>%
  as.matrix()

df_FE_X <- data %>% 
  select(state, year, lnC_it_lag, lnP_it, lnPn_it, lnY_it) %>%
  mutate(dummy_cols(data %>% select(year))) %>%
  pivot_longer(cols = -c(state, year)) %>%
  xtabs(data = ., value ~ year + state + name)

FE <- FE_own(df_FE_y, df_FE_X)
show(FE[1:4, ])

## Simulating properties of the FE estimator ----

## Simulation methods for bias correction and inference ----

#Define parameter grid
burn_t <- 25
rho_par <- c(0, 0.5, 0.9)
T_par = c(4, 10, 20, 50)
N_par = c(10, 100)

#Generate samples
generate_sample <- function(rho, N, T){ #Should be implemented with vectors
  yit <- matrix(data = 0, nrow = T+1, ncol = N)
  for (i in 1:N){
    y <- 0
      for (s in 1:burn_t){
        y <- rho*y + rnorm(1)
      }
    yit[1, i] <- y
    for (t in 1:T){
      y <- rho*y + rnorm(1)
      yit[t+1, i] <- y
    }
  }
  return(yit)
}

#Perform Monte Carlo
N_sim <- 10 #Number of simulations
result_MC <- data.frame(rho = double(), T = integer(), N = integer(), 
                        result_rho = double(), result_se = double())

for (rho in rho_par) {
  for (T in T_par) {
     for (N in N_par) {
       print(c(rho, T, N))
       rho_sim  <- c()
       se_sim   <- c()
       for (i in 1:N_sim) {
        df_MC   <- generate_MC(rho, N, T)
        df_MC_y <- df_MC[-1, ]
        df_MC_x <- array(df_MC[-(T+1), ], dim=c(T, N, 1))
        
        FE_MC   <- FE_own(df_MC_y, df_MC_x)
        rho_sim <- c(rho_sim, FE_MC[1])
        se_sim  <- c(se_sim, FE_MC[2])
       }
       result_MC <- rbind(result_MC, list(rho = rho, T = T, N = N, 
                                          result_rho = mean(rho_sim), result_se = mean(se_sim)))
    }
  }
}

result_MC$asymptotic_assumption <- result_MC$N / result_MC$T^3 
result_MC$nickel_bias           <- - (1+result_MC$rho)/result_MC$T
result_MC$nickel_bias_corrected <- result_MC$rho + result_MC$nickel_bias

### Half-panel jackknife bias-corrected FE estimator ----

### Bootstrap ----

# Part 2: Generalized Method of Moments (GMM) estimation ----
