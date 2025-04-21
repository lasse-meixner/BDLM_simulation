### RUNS A MINIMAL SIMULATION FOR TESTING with reps = 1

# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("i hope this works"))) # laura's initial seed

results <- run_simulation_parallel(
  model_type = c("BDML_b", "BDML_b2", "BDML_iw", "BDML_b2_iw", "BLRs"),
  N = 200,
  P = 100,
  setting = "noisy_fs",
  sigma = 1,
  simulation_size = 3, # NOTE: TESTING: REP = 1
  batch_size = 16,
  n_cores = 4)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
# 1.
first_plot <- get_combined_plots(results, save = TRUE)
# 2. zoomed in
zoomed_in_1 <- get_combined_plots_zoom(results, save = TRUE, zoom_in = c("BDML-Hier", "BDML-Basic", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-Split"), suffix = "_no_IW")
zoomed_in_2 <- get_combined_plots_zoom(results, save = TRUE, zoom_in = c("BDML-IW-Hier", "BDML-Hier", "BDML-IW", "BDML-Basic", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-Split"), suffix = "_IW")

