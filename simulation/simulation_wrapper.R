# source model-specific fitting functions
source("BDML_simulation_source.R")
source("BLRs_simulation_source.R")

## Main entry point for simulation study
run_simulation <- function(model_type, N, P, setting, sigma, simulation_size) {
  #' Run Simulation Study
  #'
  #' @param model_type A string specifying the model type to use. Options are "BDLM" or "BLRs".
  #'                   Details on the models can be found in the respective source files.
  #' @param N An integer specifying the number of observations.
  #' @param P An integer specifying the number of predictors.
  #' @param setting A string specifying the data generation setting.
  #'                Details on different settings can be found in the `generate_data_source.R` file.
  #' @param sigma A numeric value specifying the standard deviation of the noise.
  #' @param simulation_size An integer specifying the number of simulations to run for each setting.
  #' @return A data frame containing the combined results of all simulations.
  #'
  #' @details This function runs a simulation study for the specified model type.
  #'          It generates data based on the provided settings, fits the model,
  #'          and extracts the results. The function handles errors gracefully
  #'          and combines the results from all successful simulations into a single data frame.
  

  ## Generate simulation settings ----
  set.seed(abs(digest::digest2int("Who needs parallelization?")))
  seeds <- sample.int(.Machine$integer.max, size = simulation_size)
  sim_settings <- expand.grid(model_type = model_type, N = N, P = P, setting = setting, sigma = sigma, seed = seeds)
  
  ## Run all simulations ----
  results_list <- lapply(1:nrow(sim_settings), function(i) {
    # 1. Select model-specific functions
    if (sim_settings[i, "model_type"] == "BDML") {
      sim_iter <- sim_iter_BDML
    } else if (sim_settings[i, "model_type"] == "BLRs") {
      sim_iter <- sim_iter_BLRs
    } else {
      stop("Unknown model type: ", sim_settings[i, "model_type"])
    }
    # 2. Iteration message
    message("Running simulation ", i, " of ", nrow(sim_settings))
    # 3. Run simulation using model and data setting (i-th row)
    tryCatch({
      do.call(sim_iter, as.list(sim_settings[i, -1])) # pass all settings but model_type 
    }, error = function(e) {
      message("Error in simulation ", i, ": ", e$message)
      NULL  # Return NULL in case of error
    })
  })

  ## Filter and Combine all results ----
  results_list <- Filter(Negate(is.null), results_list)
  results <- do.call(rbind, results_list)
  
  return(results)
}