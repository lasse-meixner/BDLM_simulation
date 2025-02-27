library(tidyverse)

# for given seed, load both datasets

seed = 660232705
laura_file = file_name <- paste0("../ForLasse/sim_N200_P100sigma1_settingfixed_seed", seed, ".txt")
my_BRL = paste0("temp/sim_N200_P100sigma1_settingfixed_seed", seed, "_modelBLRs.txt")
my_BDML_b2 = paste0("temp/sim_N200_P100sigma1_settingfixed_seed", seed, "_modelBDML_b2.txt")

load_data <- function(file_path) {
  library(tidyverse)
  
  # Read in the file lines
  txt <- readLines(file_path)
  
  # Locate the section headers
  idx_X   <- grep("^X:", txt)
  idx_Y   <- grep("^Y:", txt)
  idx_A   <- grep("^A:", txt)
  idx_res <- grep("^Estimation Results:", txt)
  
  # Extract text blocks for each component (remove any blank lines)
  X_text   <- txt[(idx_X + 1):(idx_Y - 1)] %>% .[nzchar(trimws(.))]
  Y_text   <- txt[(idx_Y + 1):(idx_A - 1)] %>% .[nzchar(trimws(.))]
  A_text   <- txt[(idx_A + 1):(idx_res - 1)] %>% .[nzchar(trimws(.))]
  res_text <- txt[(idx_res + 1):length(txt)] %>% .[nzchar(trimws(.))]
  
  ### Process X
  # Remove any extraneous lines (e.g. those mentioning "omitted")
  X_text <- X_text[!grepl("omitted", X_text)]
  
  # Find header lines in the printed X block (they begin with optional whitespace and "[,")
  header_indices <- grep("^\\s*\\[,", X_text)
  
  if (length(header_indices) == 0) {
    # Assume the entire block is a single matrix printout.
    X_block <- paste(X_text, collapse = "\n")
    X <- as.matrix(read.table(text = X_block, header = TRUE, row.names = 1,
                                check.names = FALSE))
  } else {
    # If there are multiple header lines, split the block into sub-blocks.
    blocks <- list()
    for (i in seq_along(header_indices)) {
      start <- header_indices[i]
      end <- if (i < length(header_indices)) header_indices[i + 1] - 1 else length(X_text)
      block_text <- paste(X_text[start:end], collapse = "\n")
      block <- read.table(text = block_text, header = TRUE, row.names = 1,
                          check.names = FALSE)
      blocks[[i]] <- block
    }
    X <- do.call(cbind, blocks)
  }
  
  ### Process Y and A
  # The printed vectors include index markers like "[1]", so we remove them.
  Y_clean <- gsub("\\[[0-9]+\\]", "", paste(Y_text, collapse = " "))
  Y <- scan(text = Y_clean, quiet = TRUE)
  Y <- matrix(Y, ncol = 1)
  
  A_clean <- gsub("\\[[0-9]+\\]", "", paste(A_text, collapse = " "))
  A <- scan(text = A_clean, quiet = TRUE)
  A <- matrix(A, ncol = 1)
  
  ### Process Estimation Results
  # The results block may contain more than one printed table.
  # We first identify header lines (those starting with letters or underscores)
  res_headers <- grep("^[[:space:]]*[A-Za-z_]", res_text)
  
  if (length(res_headers) == 0) {
    results_df <- NULL
  } else if (length(res_headers) == 1) {
    # Only one block: read it as a single table.
    res_block <- paste(res_text, collapse = "\n")
    results_df <- read.table(text = res_block, header = TRUE, row.names = 1,
                             stringsAsFactors = FALSE)
  } else {
    # More than one block: split into blocks and then combine by columns.
    res_blocks <- list()
    for (i in seq_along(res_headers)) {
      start <- res_headers[i]
      end <- if (i < length(res_headers)) res_headers[i + 1] - 1 else length(res_text)
      block_text <- paste(res_text[start:end], collapse = "\n")
      block <- read.table(text = block_text, header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
      res_blocks[[i]] <- block
    }
    results_df <- do.call(cbind, res_blocks)
  }
  
  # Return a list containing the components
  return(list(X = X, Y = Y, A = A, results = results_df))
}


laura <- load_data(laura_file)
my_BRL <- load_data(my_BRL)
my_BDML_b2 <- load_data(my_BDML_b2)
