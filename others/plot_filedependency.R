# 1. Load required packages (install first if needed)
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("igraph",   quietly = TRUE)) install.packages("igraph")

library(stringr)
library(igraph)

# 2. List all R scripts in your target folder
files_full <- list.files(
  path       = "/Users/lauraliu/Library/CloudStorage/Dropbox/projects/sandbox/BDML/simulation",
  pattern    = "\\.[rR]$",
  recursive  = TRUE,
  full.names = TRUE
)

# Use just the basenames for graph nodes
files <- basename(files_full)

# 3. Extract all source() dependencies into a data.frame, using basenames
edges_df <- do.call(rbind, lapply(seq_along(files_full), function(i) {
  f_full <- files_full[i]
  f      <- files[i]
  lines   <- readLines(f_full, warn = FALSE)
  matches <- str_match_all(lines, "source\\s*\\(\\s*['\"]([^'\"]+)['\"]")
  srcs    <- unique(unlist(lapply(matches, function(m) m[,2]), use.names = FALSE))
  if (length(srcs) == 0) return(NULL)
  # resolve relative paths and strip to basenames
  srcs_full <- file.path(dirname(f_full), srcs)
  data.frame(
    from = f,
    to   = basename(srcs_full),
    stringsAsFactors = FALSE
  )
}))

# 4. Build an igraph object with bare filenames
g <- graph_from_data_frame(
  d        = edges_df,
  vertices = data.frame(name = files, stringsAsFactors = FALSE),
  directed = TRUE
)

# 5. Inspect and plot
print(edges_df)               
cat("Total edges:", ecount(g), "\n")
plot(
  g,
  vertex.label.cex = 0.8,
  vertex.size       = 15,
  edge.arrow.size   = 0.4
)
