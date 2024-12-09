# Loading source functions from "simulation_source.R" for data generation and simulation wrapper
source("simulation_wrapper.R")
source("simulation_wrapper_parallel.R")

# Run simulation ----
results <- run_simulation(N = 200, P = 100, model_type = c("BLRs", "BDML"), setting = c("joint"), sigma = c(1, 2, 4), simulation_size = 100)

# Note: for larger simulations, use run_simulation_parallel instead with additional arguments for batch_size and n_cores
# results <- run_simulation_parallel(N = 200, P = 100, model_type = c("BLRs", "BDML"), setting = c("joint"), sigma = c(1, 2, 4), simulation_size = 100, batch_size = 64, n_cores = 4)

## Summarize results ----
results_table <- results |>
  group_by(Method, setting, sigma, N, P) |>
  summarise(coverage = mean(catch), 
            rmse = sqrt(mean(squared_error)), 
            width = mean(interval_width))

## print results table
xtable::xtable(results_table)

# Save results to csv ----

## create results directory if it does not exist
if (!dir.exists("../results")) {
  dir.create("../results")
}

## get timestamp and save
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M") # get timestamp
results_file <- paste0("../results/results_", timestamp, ".csv")
write.csv(results_table, results_file, row.names = FALSE)

## create and save plots
plot_results(results_table, today_str, zoom = FALSE)
plot_results(results_table, today_str, zoom = TRUE)