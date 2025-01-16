# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("my simulation settings")))

results <- run_simulation_parallel(
  model_type = c("BDML_b2", "BDML_r2d2", "BLRs"),
  N = 200,
  P = 100,
  setting = "fixed",
  sigma = c(1, 2, 4),
  simulation_size = 200,
  batch_size = 4, 
  n_cores = 2) 

## Summarize results ----
results_table <- make_results_table(results) # from results_plotting_source.R

## print results table
xtable::xtable(results_table)

## create and save plots
get_combined_plots(results, save = TRUE)
get_combined_plots_zoom(results, save = TRUE)