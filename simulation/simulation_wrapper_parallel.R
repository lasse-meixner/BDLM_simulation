library(pbmcapply)
library(log4r)

# source model-specific fitting functions
source("BDML_simulation_source.R")
source("BLRs_simulation_source.R")

## Main entry point for simulation study
run_simulation_parallel <- function(model_type, n, p, R_Y2, R_D2, rho, alpha, simulation_size, 
                                    batch_size = 64, n_cores = 4, save_checkpoints = FALSE,
                                    datetime_tag = format(Sys.time(), "%Y%m%d-%H%M")) {
    #' Run a parallelized simulation study
    #'
    #' @param model_type Character vector specifying the model types to be simulated (e.g., "BDML", "BLRs").
    #' @param n Integer vector specifying the sample sizes for the simulations.
    #' @param p Integer vector specifying the number of predictors for the simulations.
    #' @param R_Y2 Numeric vector specifying the partial R2 in the Y equation.
    #' @param R_D2 Numeric vector specifying the R2 in the D equation.
    #' @param rho Numeric vector specifying the correlation between beta and gamma.
    #' @param alpha Numeric vector specifying the true value of the parameter of interest.
    #' @param simulation_size Integer specifying the total number of simulations to run per setting.
    #' @param batch_size Integer specifying the number of simulations to run in each batch. Default is 64.
    #' @param n_cores Integer specifying the number of CPU cores to use for parallel execution. Default is 4.
    #' @param save_checkpoints Logical indicating whether to save intermediate results after each batch. Default is FALSE.
    #' @param datetime_tag Character string for tagging the output files with the current date and time.
    #' @return A list of data frames, each containing the results of a batch of simulations.
    


    # Generate simulation settings
    seeds <- sample.int(.Machine$integer.max, size = simulation_size)
    sim_settings <- expand.grid(model_type = model_type, 
                                n = n, 
                                p = p, 
                                R_Y2  = R_Y2,
                                R_D2  = R_D2,
                                rho   = rho,
                                alpha = alpha, 
                                seed = seeds,
                                stringsAsFactors = FALSE) |>
                                arrange(model_type) # order by model type
    
    # init print
    print(paste0(
        "START-UP: Running a total of ", 
        (nrow(sim_settings) / length(model_type)),
        " simulations for each model, with ", 
        sum(sim_settings$model_type == "BLR-baseline") * 3, # Naive, HCHP, Linero
        " BLR-baseline (Naive, HCHP, Linero) fits, ",
        sum(sim_settings$model_type == "BLRs-FDML") * 4, # FDML-Full, FDML-Split, FDML-XFit, FDML-Alt
        " BLRs-FDML fits, ",
        sum(sim_settings$model_type == "BLRs-OLS-oracle") * 2, # OLS and oracle
        " BLRs-OLS-oracle (OLS & oracle) fits, ",
        sum(sim_settings$model_type == "BDML-IW"),
        " BDML-IW fits, ",
        sum(sim_settings$model_type == "BDML-IW-HP"),
        " BDML-IW-HP fits, ",
        sum(sim_settings$model_type == "BDML-LKJ"),
        " BDML-LKJ fits, ",
        sum(sim_settings$model_type == "BDML-LKJ-HP"), 
        " BDML-LKJ-HP fits, and ",
        sum(sim_settings$model_type == "BDML-R2D2"),
        " BDML-R2D2 fits, on ", n_cores, " cores. For details, see the log file."))

     # Initialize logger
    if (!dir.exists("logs")) dir.create("logs")
    log_file <- paste0("logs/simulation_log_", datetime_tag, ".log")
    logger <- create.logger(logfile = log_file, level = "INFO")
    
    # Log simulation start time
    start_time <- Sys.time()
    info(logger, paste("Simulation started at:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    
    # simulation settings logging
    # Log the input vectors
    info(logger, "Input vectors passed to the simulation function:")
    info(logger, paste("  model_type:", paste(model_type, collapse = ", ")))
    info(logger, paste("  n:", paste(n, collapse = ", ")))
    info(logger, paste("  p:", paste(p, collapse = ", ")))
    info(logger, paste("  R_Y2:", paste(R_Y2, collapse = ", ")))
    info(logger, paste("  R_D2:", paste(R_D2, collapse = ", ")))
    info(logger, paste("  rho:", paste(rho, collapse = ", ")))
    info(logger, paste("  alpha:", paste(alpha, collapse = ", ")))
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
            tryCatch({
                sim_iter <- switch(
                    batch_settings[i, "model_type"],
                    "BDML-LKJ" = sim_iter_bdml_lkj,
                    "BDML-LKJ-HP" = sim_iter_bdml_lkj_hp,
                    "BDML-IW" = sim_iter_bdml_iw,
                    "BDML-IW-HP" = sim_iter_bdml_iw_hp,
                    "BDML-IW-JS-MAT" = sim_iter_bdml_iw_js_mat,
                    "BDML-IW-JS-I" = sim_iter_bdml_iw_js_i,
                    "BLRs-baseline"   = sim_iter_BLRs_baseline,
                    "BLRs-FDML"       = sim_iter_BLRs_FDML,
                    "BLRs-OLS-oracle" = sim_iter_BLRs_OLS_oracle,
                    "BDML-R2D2" = sim_iter_bdml_r2d2,
                    stop("Unknown model type: ", batch_settings[i, "model_type"])
                )
            
                info(logger, paste("Running simulation", i, "in batch", batch, "with model setting", batch_settings[i, "model_type"]))
            
            # Run simulation and handle errors and warnings (Note: this was a pain since some is thrown by low-level stan execution, this was gold: https://stackoverflow.com/questions/49694552/suppress-messages-from-underlying-c-function-in-r/49722545#49722545)
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
    result_file <- paste0("results/results_", datetime_tag, ".csv")
    write.csv(results, result_file, row.names = FALSE)
    
    info(logger, paste("Saved final results to", result_file))

    # Log information about failed models
    simulations_per_model <- nrow(sim_settings) / length(model_type)
    if (length(failed_models) > 0) {
        info(logger, "FITTING SUMMARY: The following models failed:")
        for (model in names(failed_models)) {
            info(logger, paste("  ", model, ": ", failed_models[model], "failed simulations out of ", simulations_per_model, "total simulations"))
        }
    }

    # Log simulation completion time
    end_time <- Sys.time()
    info(logger, paste("Simulation completed at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "for a total duration of", round(difftime(end_time, start_time, units = "secs"), 2), "seconds"))
    
    return(results)
}
