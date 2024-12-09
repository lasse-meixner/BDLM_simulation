library(pbmcapply)

# source model-specific fitting functions
source("BDML_simulation_source.R")
source("BLRs_simulation_source.R")

## Main entry point for simulation study
run_simulation_parallel <- function(model_type, N, P, setting, sigma, simulation_size, batch_size = 64, n_cores = 4, save_checkpoints = FALSE) {
    #' Run a parallelized simulation study
    #'
    #' @param model_type Character vector specifying the model types to be simulated (e.g., "BDML", "BLRs").
    #' @param N Integer vector specifying the sample sizes for the simulations.
    #' @param P Integer vector specifying the number of predictors for the simulations.
    #' @param setting Character vector specifying different settings for the simulations.
    #' @param sigma Numeric vector specifying the noise levels for the simulations.
    #' @param simulation_size Integer specifying the total number of simulations to run.
    #' @param batch_size Integer specifying the number of simulations to run in each batch. Default is 64.
    #' @param n_cores Integer specifying the number of CPU cores to use for parallel execution. Default is 4.
    #' @return A list of data frames, each containing the results of a batch of simulations.
    

    # Generate simulation settings
    set.seed(abs(digest::digest2int("me")))
    seeds <- sample.int(.Machine$integer.max, size = simulation_size)
    sim_settings <- expand.grid(model_type = model_type, 
                                N = N, 
                                P = P, 
                                setting = setting, 
                                sigma = sigma, 
                                seed = seeds)
    
    # Determine the number of batches
    n_batches <- ceiling(nrow(sim_settings) / batch_size)
    results_list <- list()
    
    for (batch in 1:n_batches) {
        # Define the range of rows for the current batch
        start_idx <- (batch - 1) * batch_size + 1
        end_idx <- min(batch * batch_size, nrow(sim_settings))
        batch_settings <- sim_settings[start_idx:end_idx, ]
        
        # Run batch simulations in parallel
        batch_results <- pbmclapply(1:nrow(batch_settings), function(i) {
        # Select model-specific functions
        sim_iter <- switch(
            batch_settings[i, "model_type"],
            "BDML" = sim_iter_BDML,
            "BLRs" = sim_iter_BLRs,
            stop("Unknown model type: ", batch_settings[i, "model_type"])
        )
        
        # Run simulation and handle errors gracefully
        tryCatch({
            do.call(sim_iter, as.list(batch_settings[i, -1]))  # Pass all settings but model_type
        }, error = function(e) {
            message("Error in simulation ", i, ": ", e$message)
            NULL  # Return NULL on error
        })
        }, mc.cores = n_cores)  # Adjust cores for parallel execution
        
        # Filter and save intermediate results
        batch_results <- Filter(Negate(is.null), batch_results)
        results_list[[batch]] <- do.call(rbind, batch_results)
        
        # Save intermediate batch results (optional)
        if (save_checkpoints) {
            if (!dir.exists("temp")) {dir.create("temp")}
            saveRDS(results_list[[batch]], paste0("temp/batch_", batch, "_results.rds"))
        }

        message("Completed batch ", batch, " of ", n_batches)
    }
    
    # Combine all results
    results <- do.call(rbind, results_list)
    return(results)
}