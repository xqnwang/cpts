ICP.update <- function(series, nfit, ncal, alpha = 0.1,
                       w = NULL, updateAlpha = FALSE, gamma = 0.005,
                       updateMethod = c("Simple", "Momentum"), momentumBW = 0.95) {
  # Set up data
  y <- series
  T <- length(y)
  alphat <- alpha
  
  # Check the series length
  if (T < (nfit+ncal+1))
    stop("series length must be larger than nfit+ncal")
  
  # Check the weights
  if (is.null(w))
    w <- rep(1, T)
  
  # Check gamma - step size parameter
  if (is.null(gamma) || (gamma <= 0)) 
    gamma <- 0.005
  
  # Alpha update method
  updateMethod <- match.arg(updateMethod)
  
  # Initialize data storage variables
  predSeq <- loSeq <- upSeq <- rep(0, T-nfit-ncal)
  scores <- rep(0, T-nfit)
  if (updateAlpha) {
    alphaSeq <- rep(alpha, T-nfit-ncal);
    errSeq <- rep(0, T-nfit-ncal);
  }
  
  for (t in (nfit+1):T) {
    # Fit AR model and compute new conformity score
    out <- forecast::Arima(y[(t-nfit):(t-1)], order = c(2,0,0),
                           include.mean = TRUE, method = "CSS")
    fc <- forecast::forecast(out, h = 1)
    pred <- as.vector(fc$mean)
    
    # Compute conformity score
    scores[t-nfit] <- abs(y[t] - pred)
    
    # Get forecasts on test set
    if (t > (nfit+ncal)) {
      recentScores <- scores[(t-nfit-ncal):(t-nfit-1)]
      recentWeights <- w[(t-ncal):(t-1)]
      q <- weighted.quantile(c(recentScores, Inf), 
                             prob = 1-alphat, 
                             w = c(recentWeights, w[t]),
                             sorted = FALSE)
      predSeq[t-nfit-ncal] <- pred
      loSeq[t-nfit-ncal] <- pred - q
      upSeq[t-nfit-ncal] <- pred + q
      
      if (updateAlpha) {
        # Compute errt
        errSeq[t-nfit-ncal] <- as.numeric(scores[t-nfit] > quantile(recentScores, 1-alphat))
        
        # Update alphat
        alphaSeq[t-nfit-ncal] <- alphat
        if (updateMethod == "Simple") {
          alphat <- alphat + gamma*(alpha - errSeq[t-nfit-ncal])
        } else if (updateMethod == "Momentum") {
          p <- rev(momentumBW^(1:(t-nfit-ncal)))
          p <- p/sum(p)
          alphat <- alphat + gamma*(alpha - sum(errSeq[1:(t-nfit-ncal)]*p))
        }
      }
    }
  }
  
  if (!updateAlpha) {
    return(list(pred = predSeq, lo = loSeq, up = upSeq))
  } else {
    return(list(pred = predSeq, lo = loSeq, up = upSeq, alpha = alphaSeq))
  }
}
