ICP.split <- function(data, x0, alpha = 0.1, rho = 0.5, split = NULL, w = NULL) {
  
  # Set up data
  n <- nrow(data)
  p <- ncol(data)-1
  n0 <- nrow(x0)
  
  # Check the weights
  if (is.null(w)) 
    w <- rep(1, n+n0)
  
  # Check rho
  if (is.null(rho) || length(rho) != 1 || !is.numeric(rho) || rho <= 0 || rho >= 1) 
    stop("rho must be a number in between 0 and 1")
  
  # If the user passed indices for the split, use them
  if (!is.null(split)) {
    i1 <- split
  } else {
    # Otherwise make a split using rho
    i1 <- seq.int(floor(n*rho))
  }
  i2 <- (1:n)[-i1]
  n1 <- length(i1)
  n2 <- length(i2)
  
  # Initialize data storage variables
  lo <- up <- rep(0, n0)
  
  # Train on first part
  out <- lm(y~., data = data[i1, ])
  
  # Get residuals and quantiles on second
  res <- abs(as.vector(data[i2, "y"]) - as.vector(predict(out, data[i2, !(names(data) %in% c("y"))])))
  
  # Get forecasts on test set
  pred <- as.vector(predict(out, x0))
  
  for (i in 1:n0) {
    q <- weighted.quantile(c(res, Inf), prob = 1-alpha, 
                           w = c(w[i2], w[n+i]), sorted = FALSE)
    lo[i] <- pred[i] - q
    up[i] <- pred[i] + q
  }
  
  return(list(pred = pred, lo = lo, up = up, split = i1))
}

weighted.quantile <- function(v, prob, w = NULL, sorted = FALSE) {
  if (is.null(w)) w <- rep(1, length(v))
  if (!sorted) { o <- order(v); v <- v[o]; w <- w[o] }
  i <- which(cumsum(w/sum(w)) >= prob)
  if (length(i) == 0) return(Inf) # Can happen with infinite weights
  else return(v[min(i)])
}