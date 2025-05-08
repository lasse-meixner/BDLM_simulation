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

# Create and save plots
## 1. all methods
get_combined_plots(results, save = TRUE, datetime_tag = datetime_tag)
## 2. zoom in
get_combined_plots_zoom(results, save = TRUE, zoom_in = c("BDML-LKJ-HP", "BDML-LKJ", "Linero"),
                        datetime_tag = datetime_tag)
