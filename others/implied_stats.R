# Set working directory to current script location
if (rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
    # Fallback for non-RStudio environments
    setwd(dirname(sys.frame(1)$ofile))
}

library(tidyverse)
library(mvtnorm)
library(openxlsx)

# generate_data <- function(n, p, R_Y2, R_D2, rho, alpha) {

#   # 1a) build the 2×2 covariance for (beta_j,gamma_j)'
#   Sigma <- (1/p) * matrix(
#     c(R_Y2,
#       rho * sqrt(R_Y2 * R_D2),
#       rho * sqrt(R_Y2 * R_D2),
#       R_D2),
#     nrow = 2, byrow = TRUE
#   )

#   # 1b) draw p independent coefficient pairs
#   #     (β_j, γ_j)' ~ N(0, Sigma)
#   B <- rmvnorm(n = p, mean = c(0, 0), sigma = Sigma)
#   beta  <- B[,1]
#   gamma <- B[,2]

#   # 2) residual standard‐deviations
#   sigma_eps <- sqrt(1 - R_Y2)
#   sigma_V   <- sqrt(1 - R_D2)

#   # 3) simulate data
#   X   <- matrix(rnorm(n * p), nrow = n, ncol = p)        # regressors
#   V   <- rnorm(n, mean = 0, sd = sigma_V)                # first‐stage noise
#   eps <- rnorm(n, mean = 0, sd = sigma_eps)              # outcome noise

#   D <- as.numeric(X %*% gamma + V)
#   Y <- as.numeric(alpha * D + X %*% beta + eps)

#   # return a list
#   list(
#     X     = X,
#     D     = D,
#     Y     = Y,
#     R_Y2  = R_Y2,
#     R_D2  = R_D2,
#     rho   = rho,
#     alpha = alpha
#   )
# }

# Define parameter grids
R_Y2_vals <- c(0, 0.4, 0.8)
R_D2_vals <- c(0, 0.4, 0.8)
rho_vals <- c(-0.5, 0, 0.5)
alpha_vals <- c(1 / 4, 1 / 2, 1)

implied_stats <- expand_grid(
    R_Y2 = R_Y2_vals,
    R_D2 = R_D2_vals,
    rho = rho_vals,
    alpha = alpha_vals
) %>%
    mutate(
        SB = rho * sqrt(R_Y2 * R_D2) / alpha,
        ES = alpha / sqrt(1 + alpha^2 + 2 * alpha * rho * sqrt(R_Y2 * R_D2)),
        row_key = paste0(R_Y2, ",", R_D2),
        col_key = paste0(rho, ",", alpha)
    )

v_stats <- c("SB", "ES")

# Create workbook
wb <- createWorkbook()

current_row <- 1

for (v_stat in v_stats) {
    # Pivot to create 9x9 table
    implied_stats_table <- implied_stats %>%
        pivot_wider(
            names_from = col_key,
            values_from = v_stat,
            id_cols = row_key
        ) %>%
        column_to_rownames("row_key")

    # Add worksheet if it doesn't exist
    if (!"implied_stats" %in% names(wb)) {
        addWorksheet(wb, "implied_stats")
    }

    # Write header
    writeData(wb, "implied_stats",
        paste0(v_stat, " Table (R_Y2/R_D2 combinations in rows, rho/alpha combinations in columns)"),
        startCol = 1, startRow = current_row
    )
    current_row <- current_row + 1

    # Prepare table for writing
    implied_stats_table <- implied_stats_table %>%
        round(3) %>%
        rownames_to_column("R_Y2_R_D2")
    names(implied_stats_table)[1] <- ""

    # Write table
    writeData(wb, "implied_stats", implied_stats_table,
        startCol = 1, startRow = current_row, colNames = TRUE
    )
    current_row <- current_row + nrow(implied_stats_table) + 2 # +2 for header and spacing
}

# Save workbook
saveWorkbook(wb, "implied_stats.xlsx", overwrite = TRUE)
