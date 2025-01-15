# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("my simulation settings")))

results <- run_simulation_parallel(
  model_type = c("BDML", "BLRs"),
  N = 200,
  P = 100,
  setting = "fixed",
  sigma = c(1, 2, 4),
  simulation_size = 200,
  seeds = seeds
)

## Summarize results ----
results_table <- make_results_table(results)

## print results table (its saved already by the runner)
xtable::xtable(results_table)

## create and save plots
get_combined_plots(results, save = TRUE)
get_combined_plots_zoom(results, save = TRUE)