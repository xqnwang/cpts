#' Create lags of a matrix
#' 
#' Compute a lagged version of a matrix, shifting the time base back
#' by a given number of observations for each column.
#' 
#' @param x A matrix or multivariate time series.
#' @param lag A vector of lags with a length equal to that of the columns of \code{x}.
#' @return An object with the same class as \code{x}.
#' @example
#' 
#' x <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' lagmatrix(x, c(1,1,2,3))
#' 
#' @export
lagmatrix <- function(x, lag) {
  # Ensure 'x' is a matrix
  if (!is.matrix(x))
    stop("ensure x is a matrix")
  if (any(lag < 1))
    stop("ensure lag is a positive integer")
  n <- NROW(x)
  k <- length(lag)
  
  if (NCOL(x) != k)
    stop("lag must have the same number of columns as x")
  
  lmat <- x
  for (i in 1:k) {
    lmat[, i] <- c(rep(NA, lag[i]), x[1:(n - lag[i]), i])
  }
  return(structure(lmat, class = class(x)))
}
