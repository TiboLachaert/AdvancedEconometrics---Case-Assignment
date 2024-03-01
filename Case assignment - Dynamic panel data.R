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
burn_t <- 25
rho_sim <- c(0, 0.5, 0.9)
N_sim = c(10, 100)
T_sim = c(4, 10, 20, 50)

generate_MC <- function(rho, N, T){
  yit <- matrix(data = 0, nrow = N, ncol = T)
  for (i in 1:N){
    y <- 0
      for (s in 1:burn_t){
        y <- rho*y + rnorm(1)
      }
    for (t in 1:T){
      y <- rho*y + rnorm(1)
      yit[i, t] <- y
    }
  }
  return(yit)
}

for (rho in rho_sim) {
  for (N in N_sim){
    for (T in T_sim){
      print(c(rho, N, T))
    }
  }
}

### Half-panel jackknife bias-corrected FE estimator ----

### Bootstrap ----

# Part 2: Generalized Method of Moments (GMM) estimation ----
