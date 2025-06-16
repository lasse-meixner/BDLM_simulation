# Set working directory to current script location
if (rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
    # Fallback for non-RStudio environments
    setwd(dirname(sys.frame(1)$ofile))
}

library(tidyverse)
library(mvtnorm)

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
alpha_vals <- c(1 / 4, 1, 4)

AsymSB_table <- function(alpha_vals, R_D2_vals) {
    # Create combinations for AsymSB only

    asym_grid <- expand_grid(
        alpha = alpha_vals,
        R_D2 = R_D2_vals
    ) %>%
        mutate(AsymSB = alpha * (1 - R_D2))

    # Pivot to create table format
    asym_table <- asym_grid %>%
        pivot_wider(
            names_from = R_D2,
            values_from = AsymSB
        ) %>%
        column_to_rownames("alpha")

    return(asym_table)
}

# Create and display AsymSB table
out_AsymSB_table <- AsymSB_table(alpha_vals, R_D2_vals)
cat("\nAsymSB Table (alpha in rows, R_D2 in columns):\n")
print(out_AsymSB_table)

r_DY_9x9_table <- function(R_Y2_vals, R_D2_vals, rho_vals, alpha_vals) {
    # Create all combinations
    param_grid <- expand_grid(
        R_Y2 = R_Y2_vals,
        R_D2 = R_D2_vals,
        rho = rho_vals,
        alpha = alpha_vals
    )

    # Calculate r_DY
    results <- param_grid %>%
        mutate(
            r_DY = alpha^2 / (alpha^2 + 1 + 2 * alpha * rho * sqrt(R_Y2 * R_D2)),
            row_key = paste0(R_Y2, ",", R_D2),
            col_key = paste0(rho, ",", alpha)
        ) %>%
        select(row_key, col_key, r_DY)

    # Pivot to create 9x9 table
    r_DY_table <- results %>%
        pivot_wider(
            names_from = col_key,
            values_from = r_DY
        ) %>%
        column_to_rownames("row_key")

    return(r_DY_table)
}

# Create and display r_DY 9x9 table
out_r_DY_table <- r_DY_9x9_table(R_Y2_vals, R_D2_vals, rho_vals, alpha_vals)
cat("\nr_DY 9x9 Table (R_Y2/R_D2 combinations in rows, rho/alpha combinations in columns):\n")
print(round(out_r_DY_table, 3))

# Export both tables to CSV sequentially
write_lines('"AsymSB Table (alpha in rows R_D2 in columns)"', "parameter_tables.csv")
asymsb_table <- out_AsymSB_table %>% rownames_to_column("alpha")
names(asymsb_table)[1] <- ""
write_csv(asymsb_table, "parameter_tables.csv", append = TRUE, col_names = TRUE)

write_lines("", "parameter_tables.csv", append = TRUE)
write_lines('"r_DY Table (R_Y2/R_D2 combinations in rows rho/alpha combinations in columns)"', "parameter_tables.csv", append = TRUE)
rdy_table <- out_r_DY_table %>%
    round(3) %>%
    rownames_to_column("R_Y2_R_D2")
names(rdy_table)[1] <- ""
write_csv(rdy_table, "parameter_tables.csv", append = TRUE, col_names = TRUE)
