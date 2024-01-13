# The method commonly used in CP papers, which corresponds to Type 1 in Hyndman & Fan (1996).
weighted.quantile <- function(v, prob, w = NULL, sorted = FALSE) {
  if (is.null(w)) w <- rep(1, length(v))
  if (!sorted) { o <- order(v); v <- v[o]; w <- w[o] }
  i <- which(cumsum(w/sum(w)) >= prob)
  if (length(i) == 0) return(Inf) # Can happen with infinite weights
  else return(v[min(i)])
}