### RUNS A MINIMAL SIMULATION FOR TESTING with reps = 1

# Use parallelized version by default
source("simulation_wrapper_parallel.R")
source("results_plotting_source.R")

### MAIN ENTRY POINT FOR SIMULATION
# Run simulation ----
set.seed(abs(digest::digest2int("Bayesian Double Machine Learning for Causal Inference"))) 

datetime_tag <- format(Sys.time(), "%Y%m%d-%H%M")
results <- run_simulation_parallel(
  model_type = c("BLRs-baseline", "BLRs-FDML","BLRs-OLS-oracle"),
  n = 200,
  p = 100,
  R_Y2 = c(0, 0.5),
  R_D2 = c(0, 0.5),
  rho  = .5,
  alpha= 1,
  simulation_size = 2, # NOTE: low testing rep
  batch_size = 16,
  n_cores = 4,
  datetime_tag = datetime_tag)

## Summarize and print results (for curiosity - all saved to disk automatically) ----
results_table <- make_results_table(results) # from results_plotting_source.R
xtable::xtable(results_table)

## create and save plots
grid <- expand.grid(
  R_D2 = unique(results$R_D2),
  rho  = unique(results$rho),
  alpha= unique(results$alpha)
)

for(i in seq_len(nrow(grid))){
  args <- grid[i, ]
  subres <- results %>%
    filter(R_D2 == args$R_D2,
           rho   == args$rho,
           alpha == args$alpha)
  
  spec_tag = paste0(datetime_tag, "_RD2_", args$R_D2, "_rho_", args$rho,
                          "_alpha_", args$alpha)
  
  # 1. all
  first_plot <- get_combined_plots(results = subres, save = file.path("results", datetime_tag), suffix = spec_tag)
  # 2. zoomed in
  zoomed_in <- get_combined_plots_zoom(results = subres, save = file.path("results", datetime_tag), 
                                         zoom_in = c("BDML-LKJ-HP", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-XFit", "FDML-Alt", "OLS", "Oracle"),
                                         suffix = paste0(spec_tag, "_top_coverage"))
  }