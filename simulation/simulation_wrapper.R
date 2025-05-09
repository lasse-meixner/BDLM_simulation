# source model-specific fitting functions
source("BDML_simulation_source.R")
source("BLRs_simulation_source.R")

## Main entry point for simulation study
run_simulation <- function(model_type, n, p, R_Y2, R_D2, rho, alpha, simulation_size, datetime_tag) {
  #' Run Simulation Study
  #'
  #' @param model_type A string specifying the model type to use. Options are "BDML_b2", "BDML_r2d2", and "BLRs".
  #'                   Details on the models can be found in the respective source files.
  #' @param n An integer specifying the number of observations.
  #' @param p An integer specifying the number of predictors.
  #' @param R_Y2 Numeric vector specifying the partial R2 in the Y equation.
  #' @param R_D2 Numeric vector specifying the R2 in the D equation.
  #' @param rho Numeric vector specifying the correlation between beta and gamma.
  #' @param alpha Numeric vector specifying the true value of the parameter of interest.
  #' @param simulation_size An integer specifying the number of simulations to run for each setting.
  #' @param datetime_tag Character string for tagging the output files with the current date and time.
  #' @return A data frame containing the combined results of all simulations.
  #'
  #' @details This function runs a simulation study for the specified model type.
  #'          It generates data based on the provided settings, fits the model,
  #'          and extracts the results. The function handles errors gracefully
  #'          and combines the results from all successful simulations into a single data frame.
  

  ## Generate simulation settings ----
  set.seed(abs(digest::digest2int("Bayesian Double Machine Learning for Causal Inference")))
  seeds <- sample.int(.Machine$integer.max, size = simulation_size)
  sim_settings <- expand.grid(model_type = model_type, n = n, p = p, R_Y2 = R_Y2, R_D2 = R_D2, rho = rho, alpha = alpha, seed = seeds)
  
  ## Run all simulations ----
  results_list <- lapply(1:nrow(sim_settings), function(i) {
    # 1. Select model-specific functions
    if (sim_settings[i, "model_type"] == "BDML-LKJ-HP") {
      sim_iter <- sim_iter_bdml_lkj_hp
    } else if (sim_settings[i, "model_type"] == "BDML-R2D2") {
      sim_iter <- sim_iter_bdml_r2d2
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

  ## save results to csv
  write.csv(results, paste0("results/results_", datetime_tag, ".csv"), row.names = FALSE)
  
  return(results)
}
