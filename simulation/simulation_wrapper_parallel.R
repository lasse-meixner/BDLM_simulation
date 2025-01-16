library(pbmcapply)
library(log4r)

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
    #' @param simulation_size Integer specifying the total number of simulations to run per setting.
    #' @param batch_size Integer specifying the number of simulations to run in each batch. Default is 64.
    #' @param n_cores Integer specifying the number of CPU cores to use for parallel execution. Default is 4.
    #' @return A list of data frames, each containing the results of a batch of simulations.
    

     # Initialize logger
    if (!dir.exists("logs")) dir.create("logs")
    log_file <- paste0("logs/simulation_log_", format(Sys.time(), "%Y%m%d%H%M%S"), ".log")
    logger <- create.logger(logfile = log_file, level = "INFO")
    
    info(logger, paste("Simulation started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    # Generate simulation settings
    set.seed(abs(digest::digest2int("me")))
    seeds <- sample.int(.Machine$integer.max, size = simulation_size)
    sim_settings <- expand.grid(model_type = model_type, 
                                N = N, 
                                P = P, 
                                setting = setting, 
                                sigma = sigma, 
                                seed = seeds,
                                stringsAsFactors = FALSE)

    # simulation settings logging
    # Log the input vectors
    info(logger, "Input vectors passed to the simulation function:")
    info(logger, paste("  model_type:", paste(model_type, collapse = ", ")))
    info(logger, paste("  N:", paste(N, collapse = ", ")))
    info(logger, paste("  P:", paste(P, collapse = ", ")))
    info(logger, paste("  setting:", paste(setting, collapse = ", ")))
    info(logger, paste("  sigma:", paste(sigma, collapse = ", ")))
    info(logger, paste("  simulation_size:", simulation_size))
    info(logger, paste("  batch_size:", batch_size))
    info(logger, paste("  n_cores:", n_cores))
    info(logger, paste("  save_checkpoints:", ifelse(save_checkpoints, "TRUE", "FALSE")))
    info(logger, paste("Generated", nrow(sim_settings), "simulation rows."))

    # init print
    print(paste0(
        "Running a total of ", 
        (nrow(sim_settings) / length(model_type)),
        " simulations for each model, with # BLRs: ", 
        sum(sim_settings$model_type == "BLRs")*5, 
        " and # BDML: ", 
        sum(sim_settings$model_type == "BDML")*2, 
        " on ", n_cores, " cores.\n"))
    
    # Determine the number of batches
    n_batches <- ceiling(nrow(sim_settings) / batch_size)
    results_list <- list()
    
    for (batch in 1:n_batches) {
        batch_start_time <- Sys.time()
        info(logger, paste("Starting batch", batch, "of", n_batches))
        
        # Define the range of rows for the current batch
        start_idx <- (batch - 1) * batch_size + 1
        end_idx <- min(batch * batch_size, nrow(sim_settings))
        batch_settings <- sim_settings[start_idx:end_idx, ]
        
        # Run batch simulations in parallel
        batch_results <- pbmclapply(1:nrow(batch_settings), function(i) {
            sim_iter <- switch(
                batch_settings[i, "model_type"],
                "BDML" = sim_iter_BDML,
                "BLRs" = sim_iter_BLRs,
                stop("Unknown model type: ", batch_settings[i, "model_type"])
            )
            
        # Run simulation and handle errors and warnings
        tryCatch({
            # Use withCallingHandlers to log warnings
            result <- withCallingHandlers({
                do.call(sim_iter, as.list(batch_settings[i, -1]))  # pass settings to model-family-specific simulator function (excluding model_type)
            }, warning = function(w) {
                # Log warnings
                warning_message <- paste("Warning in simulation", i, "in batch", batch, ":", conditionMessage(w))
                warn(logger, warning_message)
                invokeRestart("muffleWarning")  # Suppress further propagation
            })
        
            # Return result
            result
            }, error = function(e) {
                error_message <- paste("Error in simulation", i, "in batch", batch, ":", e$message)
                log4r::error(logger, error_message) 
                NULL
            })
        }, mc.cores = n_cores)  # Adjust cores for parallel execution
        
        # Filter and save intermediate results
        batch_results <- Filter(Negate(is.null), batch_results)
        results_list[[batch]] <- do.call(rbind, batch_results)
        
        # Log batch completion time
        batch_end_time <- Sys.time()
        info(logger, paste("Completed batch", batch, "in", round(difftime(batch_end_time, batch_start_time, units = "secs"), 2), "seconds"))

        
        # Save intermediate batch results (optional)
        if (save_checkpoints) {
            if (!dir.exists("temp")) {dir.create("temp")}
            saveRDS(results_list[[batch]], paste0("temp/batch_", batch, "_results.rds"))
            info(logger, "Saved checkpoint for batch %d", batch)
        }
    }
    
    # Combine all results
    results <- do.call(rbind, results_list)
    
    # Save results to CSV
    result_file <- paste0("results/results_", format(Sys.time(), "%Y%m%d%H%M"), ".csv")
    write.csv(results, result_file, row.names = FALSE)
    
    info(logger, paste("Saved final results to", result_file))
    info(logger, paste("Simulation completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


    return(results)
}