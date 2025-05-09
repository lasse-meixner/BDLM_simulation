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
  R_Y2 = c(0, 0.4, 0.8),
  R_D2 = c(0, 0.4, 0.8),
  rho  = c(-0.5, 0, 0.5),
  alpha= c(0.25, 1, 4),
  simulation_size = 2000,
  batch_size = 48,
  n_cores = 24,
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
  
  datetimespec_tag = paste0(datatime_tag, "_RD2_", args$R_D2, "_rho_", args$rho,
                          "_alpha_", args$alpha)
  # 1.
  first_plot <- get_combined_plots(results = subres, save = TRUE, datetimespec_tag = datetimespec_tag)
  # 2. zoomed in
  zoomed_in_1 <- get_combined_plots_zoom(results = subres, save = TRUE, 
                                         zoom_in = c("BDML-LKJ-HP", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-XFit", "OLS"),
                                         datetimespec_tag = datetimespec_tag)
  zoomed_in_2 <- get_combined_plots_zoom(results = subres, save = TRUE,
                                         zoom_in = c("BDML-IW-HP", "BDML-LKJ-HP", "BDML-IW", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-XFit", "OLS"),
                                         datetimespec_tag = datetimespec_tag)
}

