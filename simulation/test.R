### RUNS A MINIMAL SIMULATION FOR TESTING with reps = 1

# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("i hope this works"))) # laura's initial seed

results <- run_simulation_parallel(
  model_type = c("BDML_iw", "BDML_iw_js", "BLRs"),
  N = 200,
  P = 100,
  setting = "fixed",
  sigma = c(1, 2, 4),
  simulation_size = 10,
  batch_size = 48,
  n_cores = 24)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
# 1.
first_plot <- get_combined_plots(results, save = TRUE)
# 2. zoomed in
zoomed_in_1 <- get_combined_plots_zoom(results, save = TRUE, zoom_in = c("BDML-IW-JS-Mat", "BDML-IW-JS-I", "BDML-IW", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-Split"))

