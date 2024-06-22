# Load necessary libraries
library(plm)

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

# Function to split panel data into two halves
create_half_panel <- function(panel_data) {
  # Convert time column to numeric if it is a factor
  panel_data$time <- as.numeric(as.character(panel_data$time))
  
  # Find the maximum time value
  max_time <- max(panel_data$time)
  
  # Check if the maximum time value is even
  if (max_time %% 2 != 0) {
    stop("The time column does not have an even number of values.")
  }
  
  # Calculate the midpoint to split the data
  midpoint <- max_time / 2
  
  # Split the data into two halves
  first_half <- panel_data[panel_data$time <= midpoint, ]
  second_half <- panel_data[panel_data$time > midpoint, ]
  
  return(list(first_half = first_half, second_half = second_half))
}

# Generate panel data with N=100, T=20
rho <- 0.8
N <- 100
T <- 20
panel_data <- generate_ar1_panel(rho, N, T)

# Split the panel data into two halves
split_panels <- create_half_panel(panel_data)
first_half_panel <- split_panels$first_half
second_half_panel <- split_panels$second_half

# Function to estimate fixed effects
estimate_fixed_effects <- function(panel_data) {
  # Convert individual and time to factors
  panel_data$individual <- as.factor(panel_data$individual)
  panel_data$time <- as.factor(panel_data$time)
  
  # Fit the fixed effects model
  model <- plm(value ~ lag(value, 1), data = panel_data, 
               index = c("individual", "time"), model = "within")
  
  # Extract coefficient and Std. Error
  sum_model = summary(model)
  return(c(sum_model$coefficients[1], sum_model$coefficients[2]))
}



# Estimate fixed effects for each half-panel
fe_first_half <- estimate_fixed_effects(first_half_panel)[1]
fe_second_half <- estimate_fixed_effects(second_half_panel)[1]

# Compute the bias-corrected estimate using half-panel jackknife
fe_full <- estimate_fixed_effects(panel_data)[1]
fe_jackknife <- 2 * fe_full - 0.5 * (fe_first_half + fe_second_half)

# Output the results
cat("Fixed Effects Estimate (Full Panel):", fe_full, "\n")
cat("Fixed Effects Estimate (First Half):", fe_first_half, "\n")
cat("Fixed Effects Estimate (Second Half):", fe_second_half, "\n")
cat("Bias-Corrected Estimate (Jackknife):", fe_jackknife, "\n")





####################
# SET THE SEED
####################
set.seed(432432)


####### BLOCK BOOSTRAP WITH RENAMING THE INDIVIDUALS COLUMNS
# Block bootstrap function with new unique IDs
block_bootstrap <- function(individual_col, df) {
  # Obtain unique individuals from dataframe column
  unique_ind <- unique(individual_col)
  
  # Sample each individual to maintain persistence in the data (with replacement)
  sample_ind <- sample(unique_ind, size = length(unique_ind), replace = TRUE)
  
  # Generate new panel from this cross-section sampling
  out_df <- do.call(rbind, lapply(sample_ind, function(x) df[df$individual == x, ]))
  
  # Assign new unique IDs to each sampled individual
  out_df$individual <- as.factor(rep(1:length(sample_ind), each = nrow(df) / length(unique_ind)))
  
  return(out_df)
}

#panel_data <- generate_ar1_panel(rho=0.5, N=5, T=1)
#boot_df = block_bootstrap(panel_data$individual, panel_data)










# Perform the procedure for multiple values of N and T
rho_values <- c(0, 0.5, 0.9)
N_values <- c(10, 100)
T_values <- c(4, 10, 20, 50)
B_values <- 500 # No. of re-samples drawn

# Set value of rho
#rho = 0.8

result_MC <- data.frame(rho = double(), T = integer(), N = integer(), 
                        result_rho = double(), result_se = double(), result_se_FE = double())

for (rho in rho_values) {
  cat("RHO:", rho, "\n")
  for (N in N_values){
    cat("N:", N, "\n")
    for (T in T_values){
      cat("T:", T, "\n")
      # Generate panel data
      panel_data <- generate_ar1_panel(rho, N, T)
      estimator_list = c()
      error_list_FE = c()
      
      B <- 1
      while (B <= B_values) {
        
        #cat("B", B)
        
        panel_data_BOOT <- block_bootstrap(panel_data$individual, panel_data)
        
        # Split the panel data into two halves
        split_panels <- create_half_panel(panel_data_BOOT)
        first_half_panel <- split_panels$first_half
        second_half_panel <- split_panels$second_half
        
        # Estimate fixed effects for each half-panel
        fe_first_half <- estimate_fixed_effects(first_half_panel)[1]
        fe_second_half <- estimate_fixed_effects(second_half_panel)[1]
        
        # Compute the bias-corrected estimate using half-panel jackknife
        fe_full <- estimate_fixed_effects(panel_data)
        fe_jackknife <- 2 * fe_full[1] - 0.5 * (fe_first_half + fe_second_half)
        
        # Attach the jackknifed estimator
        #append(estimator_list, fe_jackknife)
        estimator_list = c(estimator_list, fe_jackknife)
        error_list_FE = c(fe_full[2])
        
        B <- B+1
      }
      
      result_MC <- rbind(result_MC, list(rho = rho, T = T, N = N, 
                                         result_rho = mean(estimator_list), 
                                         result_se = sd(estimator_list),    #se^b(\hat{rho}^b_FE)
                                         result_se_FE=mean(error_list_FE))) #\hat{se}(\hat{rho}^b_FE)
      
    }
  }
}
















