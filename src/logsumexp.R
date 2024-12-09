logsumexp <- function(x) {
  M <- max(x)
  return(M + log(sum(exp(x - M))))
}