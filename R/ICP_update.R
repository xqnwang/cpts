ICP.update <- function(series, nfit, ncal, alpha = 0.1,
                       weight = c("Equal", "Exp", "GLM", "RF"), base = 0.99,
                       updateAlpha = FALSE, gamma = 0.005,
                       updateMethod = c("Simple", "Momentum"), momentumBW = 0.95,
                       kess = FALSE, type = 1) {
  # Set up data
  y <- series
  T <- length(y)
  alphat <- alpha
  
  # Check the series length
  if (T < (nfit+ncal+1))
    stop("series length must be larger than nfit+ncal")
  
  # Check the weights
  weight <- match.arg(weight)
  
  # Check gamma - step size parameter
  if (updateAlpha && (is.null(gamma) || (gamma <= 0)))
    gamma <- 0.005
  
  # Alpha update method
  updateMethod <- match.arg(updateMethod)
  
  # Kish's effective sample size
  if (kess) {
    kess <- function(w) sum(w)^2 / sum(w^2)
  } else {
    kess <- NULL
  }
  
  # Initialize data storage variables
  predSeq <- loSeq <- upSeq <- rep(0, T-nfit-ncal)
  scores <- rep(0, T-nfit)
  if (updateAlpha) {
    alphaSeq <- rep(alpha, T-nfit-ncal);
    errSeq <- rep(0, T-nfit-ncal);
  }
  
  for (t in (nfit+1):T) {
    # Fit AR model and get point forecasts
    out <- forecast::Arima(y[(t-nfit):(t-1)], order = c(2,0,0),
                           include.mean = TRUE, method = "CSS")
    fc <- forecast::forecast(out, h = 1)
    pred <- as.vector(fc$mean)
    
    # Compute conformity score
    scores[t-nfit] <- abs(y[t] - pred)
    
    # Implement conformal prediction on test set
    if (t > (nfit+ncal)) {
      if (weight == "Equal") {
        recentWeights <- rep(1, ncal+1)
      } else if (weight == "Exp") {
        if (is.null(base)) base <- 0.99
        recentWeights <- c(0.99^(ncal+1-((1:ncal))), 1)
      } else if (weight %in% c("GLM", "RF")) {
        foo <- function(k) {
          partseries <- series[(t-ncal-nfit+1):t]
          c(partseries[k:(k+ncal)])
        }
        zy <- as.factor(c(rep(0, ncal), 1))
        zx <- sapply(1:nfit, foo) |> as.matrix()
        
        if (weight == "GLM") {
          obj.glm <- glm(zy ~ zx, family = "binomial")
          prob.glm <- predict(obj.glm, type = "response")
          recentWeights <- prob.glm / (1-prob.glm)
        } else if (weight == "RF") {
          obj.rf <- randomForest::randomForest(zx, zy)
          prob.rf <- predict(obj.rf, type = "prob")[,2]
          prob.rf <- pmax(pmin(prob.rf, 0.99), 0.01)
          recentWeights <- prob.rf / (1-prob.rf)
        }
      }
      
      recentScores <- scores[(t-nfit-ncal):(t-nfit-1)]
      q <- ggdist::weighted_quantile(x = c(recentScores, Inf), probs = 1-alphat,
                                     weights = recentWeights,
                                     n = kess, type = type)
      # q <- weighted.quantile(c(recentScores, Inf),
      #                        prob = 1-alphat,
      #                        w = recentWeights,
      #                        sorted = FALSE)
      predSeq[t-nfit-ncal] <- pred
      loSeq[t-nfit-ncal] <- pred - q
      upSeq[t-nfit-ncal] <- pred + q
      
      if (updateAlpha) {
        # Compute errt
        if (alphat >= 1) {
          errSeq[t-nfit-ncal] <- 1
        } else if (alphat <= 0) {
          errSeq[t-nfit-ncal] <- 0
        } else {
          errSeq[t-nfit-ncal] <- as.numeric(scores[t-nfit] > quantile(recentScores, 1-alphat))
        }
        
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

AR2.update <- function(series, nfit, ncal, alpha = 0.1) {
  y <- series
  T <- length(y)
  
  predSeq <- loSeq <- upSeq <- rep(0, T-nfit-ncal)
  for (t in (nfit+ncal+1):length(series)) {
    out <- forecast::Arima(y[(t-nfit):(t-1)], order = c(2,0,0),
                           include.mean = TRUE, method = "CSS")
    fc <- forecast::forecast(out, h = 1, level = 1-alpha)
    predSeq[t-nfit-ncal] <- as.vector(fc$mean)
    loSeq[t-nfit-ncal] <- as.vector(fc$lower)
    upSeq[t-nfit-ncal] <- as.vector(fc$upper)
  }
  
  return(list(pred = predSeq, lo = loSeq, up = upSeq))
}
