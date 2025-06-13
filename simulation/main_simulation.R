# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("Bayesian Double Machine Learning for Causal Inference"))) 

datetime_tag <- format(Sys.time(), "%Y%m%d-%H%M")

results <- run_simulation_parallel(
  model_type = c("BLRs-baseline", "BLRs-FDML","BLRs-OLS-oracle"), # available list: "BDML-LKJ","BDML-LKJ-HP","BDML-IW","BDML-IW-HP","BDML-IW-JS-MAT","BDML-IW-JS-I","BLRs-baseline","BLRs-FDML","BLRs-OLS-oracle"
  n = 200,
  p = 100,
  R_Y2 = c(0, 0.4, 0.8),
  R_D2 = 0,
  rho  = c(-0.5, 0, 0.5),
  alpha= c(0.25, 1, 4),
  simulation_size = 500,
  batch_size = 48,
  n_cores = 24,
  datetime_tag = datetime_tag)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
# 1. all
generate_stacked_plots(results, datetime_tag)
# 2. zoomed in
generate_zoomed_in_plots(results, datetime_tag, zoom_in = c("BDML-IW-JS-MAT", "BDML-IW-JS-I", "BDML-IW", "Linero"))

