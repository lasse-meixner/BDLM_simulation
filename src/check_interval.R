check_interval <- function(x, interval) {
  if(x > interval[1]) {
    if(x < interval[2])
      return(TRUE)
  }
  return(FALSE)
}