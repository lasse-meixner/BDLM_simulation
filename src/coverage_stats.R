coverage_stats <- function(theta, theta_hat, q) {
  quants         <- quantile(theta_hat, q)
  catch          <- (theta < quants[2]) & (theta > quants[1])
  interval_width <- quants[2] - quants[1]
  standard_error <- sd(theta_hat)
  point_est      <- mean(theta_hat)
  error          <- point_est - theta
  
  out <- data.frame(
    catch = catch,
    width = interval_width,
    standard_error = standard_error,
    point_est = point_est,
    error = error,
    theta = theta
  )
  rownames(out) <- NULL
  return(out)
}