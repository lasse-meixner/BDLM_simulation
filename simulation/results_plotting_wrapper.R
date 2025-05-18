source("results_plotting_source.R")

### Main entry point for results plotting [takes most recent results file and creates and saves plots]

# Load results data to plot
## Find the most recent results file path
detect_most_recent_datetime <- function(directory = "results") {
  # List all files in the specified directory
  files <- list.files(directory, pattern = "^results_\\d{12}\\.csv$", full.names = FALSE)
  
  # Check if there are any matching files
  if (length(files) == 0) {
    stop("No results files found in the directory.")
  }
  
  # Extract timestamps from filenames
  timestamps <- sub("^results_(\\d{12})\\.csv$", "\\1", files)
  
  # Convert timestamps to datetime objects for comparison
  datetimes <- as.POSIXct(timestamps, format = "%Y%m%d-%H", tz = "UTC")
  
  # Find the filename with the most recent timestamp
  most_recent_datetime <- timestamps[which.max(datetimes)]
  
  return(most_recent_datetime)
}

## load the most recent results
datetime_tag <- detect_most_recent_datetime()
results <- read.csv(paste0("results/results_", datetime_tag,".csv"), stringsAsFactors = FALSE)

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
  first_plot <- get_combined_plots(results = subres, save = file.path("results", datetime_tag), datetimespec_tag = datetimespec_tag)
  # 2. zoomed in
  zoomed_in_1 <- get_combined_plots_zoom(results = subres, save = file.path("results", datetime_tag), 
                                         zoom_in = c("BDML-LKJ-HP", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-XFit", "FDML-Alt", "OLS", "Oracle"),
                                         datetimespec_tag = datetimespec_tag)
  zoomed_in_2 <- get_combined_plots_zoom(results = subres, save = file.path("results", datetime_tag),
                                         zoom_in = c("BDML-IW-HP", "BDML-LKJ-HP", "BDML-IW", "BDML-LKJ", "Linero", "HCPH", "Naive", "FDML-Full", "FDML-XFit", "FDML-Alt", "OLS", "Oracle"),
                                         datetimespec_tag = datetimespec_tag)
}



