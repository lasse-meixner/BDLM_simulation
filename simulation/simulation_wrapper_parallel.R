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
    


    # Generate simulation settings
    seeds <- sample.int(.Machine$integer.max, size = simulation_size)
    sim_settings <- expand.grid(model_type = model_type, 
                                N = N, 
                                P = P, 
                                setting = setting, 
                                sigma = sigma, 
                                seed = seeds,
                                stringsAsFactors = FALSE) |>
                                arrange(model_type) # order by model type
    
    # init print
    print(paste0(
        "START-UP: Running a total of ", 
        (nrow(sim_settings) / length(model_type)),
        " simulations for each model, with ", 
        sum(sim_settings$model_type == "BLRs")*5, 
        " BLR fits, ",
        sum(sim_settings$model_type == "BDML_iw"),
        " BDML_iw fits, ",
        sum(sim_settings$model_type == "BDML_b2_iw"),
        " BDML_b2_iw fits, ",
        sum(sim_settings$model_type == "BDML_b"),
        " BDML_b fits, ",
        sum(sim_settings$model_type == "BDML_b2"), 
        " BDML_b2 fits, and ",
        sum(sim_settings$model_type == "BDML_r2d2"),
        " BDML_r2d2 fits, on ", n_cores, " cores. For details, see the log file."))

     # Initialize logger
    if (!dir.exists("logs")) dir.create("logs")
    log_file <- paste0("logs/simulation_log_", format(Sys.time(), "%Y%m%d%H%M%S"), ".log")
    logger <- create.logger(logfile = log_file, level = "INFO")
    
    # Log simulation start time
    start_time <- Sys.time()
    info(logger, paste("Simulation started at:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    
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

    # keep track of failed model types
    failed_models <- setNames(rep(0, length(model_type)), model_type)
    
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
            # select the appropriate simulation function based on model type
            sim_iter <- switch(
                batch_settings[i, "model_type"],
                "BDML_b" = sim_iter_BDML_b,
                "BDML_b2" = sim_iter_BDML_b2,
                "BDML_r2d2" = sim_iter_BDML_r2d2,
                "BDML_iw" = sim_iter_BDML_iw,
                "BDML_b2_iw" = sim_iter_BDML_b2_iw,
                "BLRs" = sim_iter_BLRs,
                stop("Unknown model type: ", batch_settings[i, "model_type"])
            )
        
        info(logger, paste("Running simulation", i, "in batch", batch, "with model setting", batch_settings[i, "model_type"]))
        
        # Run simulation and handle errors and warnings (Note: this was a pain since some is thrown by low-level stan execution, this was gold: https://stackoverflow.com/questions/49694552/suppress-messages-from-underlying-c-function-in-r/49722545#49722545)
        tryCatch({
            # compute simulation results
            output_str <- capture.output({
                result <- do.call(sim_iter, as.list(batch_settings[i, -1]))  # pass settings to model-family-specific simulator function (excluding model_type)
            }, type = "message")
            # log warning from stderr (lower level STAN C++ warnings)
            if (length((output_str)) > 0) {
                lapply(output_str, function(line) warn(logger, paste("STDERR output during simulation", i, 
                                                                  "in batch", batch, ":", 
                                                                  "from model", batch_settings[i, "model_type"], ":", line)))
            }
            # return result
            list(error = FALSE, result = result, model_type = batch_settings[i, "model_type"])
            }, warning = function(w) {
                # Log warnings
                warning_message <- paste("Warning in simulation", i, "in batch", batch, ":", "from model", batch_settings[i, "model_type"], ":", w$message)
                log4r::warn(logger, warning_message)
                # return result
                list(error = FALSE, result = result, model_type = batch_settings[i, "model_type"])
            }, error = function(e) {
                # Log errors
                error_message <- paste("Error in simulation", i, "in batch", batch, ":", "from model", batch_settings[i, "model_type"], ":", e$message)
                log4r::error(logger, error_message)
                # return NULL
                list(error = TRUE, result = NULL, model_type = batch_settings[i, "model_type"])
            })
        }, mc.cores = n_cores)  # Adjust cores for parallel execution
        
        # Check for errors and update failed_models
        for (j in 1:nrow(batch_settings)) {
            if (batch_results[[j]]$error) {
                failed_models[[batch_results[[j]]$model_type]] <- failed_models[[batch_results[[j]]$model_type]] + 1
            }
        }

        # Filter and save intermediate results
        batch_results <- Filter(Negate(is.null), lapply(batch_results, function(x) x$result))
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
    if (!dir.exists("results")) {dir.create("results")}
    result_file <- paste0("results/results_", format(Sys.time(), "%Y%m%d-%H%M"), ".csv")
    write.csv(results, result_file, row.names = FALSE)
    
    info(logger, paste("Saved final results to", result_file))

    # Log information about failed models
    if (length(failed_models) > 0) {
        info(logger, "FITTING SUMMARY: The following models failed:")
        for (model in names(failed_models)) {
            info(logger, paste("  ", model, ": ", failed_models[model], "failed simulations out of ", simulation_size, "total simulations"))
        }
    }

    # Log simulation completion time
    end_time <- Sys.time()
    info(logger, paste("Simulation completed at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "for a total duration of", round(difftime(end_time, start_time, units = "secs"), 2), "seconds"))
    
    return(results)
}