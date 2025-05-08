### RUNS A MINIMAL SIMULATION FOR TESTING with reps = 1

# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("Bayesian Double Machine Learning for Causal Inference"))) 

datetime_tag <- format(Sys.time(), "%Y%m%d-%H%M")
results <- run_simulation_parallel(
  model_type = c("BDML-LKJ", "BDML-LKJ-HP", "BDML-IW", "BDML-IW-HP", "BLRs"),
  n = 200,
  p = 100,
  setting = "noisy_fs",
  sigma = 1,
  simulation_size = 3, # NOTE: TESTING: REP = 1
  batch_size = 16,
  n_cores = 4,
  datetime_tag = datetime_tag)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
# 1.
first_plot <- get_combined_plots(results, save = TRUE, datetime_tag = datetime_tag)
# 2. zoomed in
zoomed_in_1 <- get_combined_plots_zoom(results, save = TRUE,
                                       zoom_in = c("BDML-LKJ-HP", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-Split"),
                                       datetime_tag = datetime_tag)
zoomed_in_2 <- get_combined_plots_zoom(results, save = TRUE, 
                                       zoom_in = c("BDML-IW-HP", "BDML-LKJ-HP", "BDML-IW", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-Split"),
                                       datetime_tag = datetime_tag)
