### ENTRY POINT TO RUN BDML VS LINERO

# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("Bayesian Double Machine Learning for Causal Inference")))

datetime_tag <- format(Sys.time(), "%Y%m%d-%H%M")
results <- run_simulation_parallel(
  # Available: "BDML-LKJ", "BDML-LKJ-HP", "BDML-IW", "BDML-IW-HP", "BDML-IW-JS-MAT", "BDML-IW-JS-I", "BLRs-baseline", "BLRs-FDML", "BLRs-OLS-oracle"
  model_type = c("BDML-IW-JS-MAT", "BLRs-baseline"),
  n = 200,
  p = 100,
  R_Y2 = c(0.05, 0.5, 0.95),
  R_D2 = c(0.05, 0.5, 0.95),
  rho = c(0.05, 0.5, 0.95),
  alpha = 0.25,
  simulation_size = 100,
  batch_size = 100,
  n_cores = 8,
  datetime_tag = datetime_tag
)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
# 1. all
generate_stacked_plots(results, datetime_tag)
# 2. zoomed in
generate_stacked_plots(results, datetime_tag, zoom_in = c("BDML-IW-JS-MAT", "BLRs-baseline"))

## RMSE tablelibrary(dplyr)
library(tidyr)

# Define the parameters to match on (adjust as needed)
params <- c("rho", "alpha", "n", "p", "R_Y2", "R_D2")

# Calculate RMSE for each method and parameter setting
rmse_df <- results %>%
  group_by(across(all_of(params)), Method) %>%
  summarise(RMSE = sqrt(mean(squared_error, na.rm = TRUE)), .groups = "drop")

# Reshape to wide format for comparison
rmse_wide <- rmse_df %>%
  pivot_wider(names_from = Method, values_from = RMSE)

# Find cases where linero's RMSE is larger than BDML's
print(rmse_wide, n = Inf)
