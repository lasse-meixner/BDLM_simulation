autocorr.mat <- function(p = 100, rho = 0.99) {
  mat <- diag(p)
  autocorr <- sapply(seq_len(p), function(i) rho^i)
  mat[lower.tri(mat)] <- autocorr[unlist(sapply(seq.int(p-1, 1), function(i) seq_len(i)))]
  mat[upper.tri(mat)] <- autocorr[unlist(sapply(seq_len(p - 1), function(i) seq.int(i, 1)))]
  return(mat)
}