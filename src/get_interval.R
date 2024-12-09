get_interval <- function(x) {
  mu <- x$bF[1]
  sigma <- x$SD.bF[1]
  return(mu + c(-1.96,1.96) * sigma)
}