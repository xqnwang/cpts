PID.update <- function(series, nfit, nburnin = 100, alpha = 0.1,
                       integrate = TRUE, scorecast = TRUE,
                       lr = 0.1, Csat = 1, KI = 200) {
  # Set up data
  y <- series
  T <- length(y)
  
  # Check the series length
  if (T < (nfit+nburnin+1))
    stop("series length must be larger than nfit+nburnin")
  
  # Initialize data storage variables
  predSeq <- loSeq <- upSeq <- rep(0, T-nfit)
  scores <- errs <- qs <- qts <- integrators <- scorecasts <- rep(0, T-nfit)
  
  for (t in (nfit+1):T) {
    # Learning rate
    t_lr <- t-fit-1
    t_lr_min <- max(t_lr-nburnin+1, 1)
    lr_t <- ifelse(t_lr > 0, lr*(max(scores[t_lr_min, t_lr]) - min(scores[t_lr_min, t_lr])), lr)
    
    # Fit AR model and get point forecasts
    out <- forecast::Arima(y[(t-nfit):(t-1)], order = c(2,0,0),
                           include.mean = TRUE, method = "CSS")
    fc <- forecast::forecast(out, h = 1)
    pred <- as.vector(fc$mean)
    
    # Observe y_t, calculate conformity score and err
    scores[t-nfit] <- abs(y[t] - pred)
    errs[t-nfit] <- !(scores[t-nfit] <= qs[t-nfit])
    
    # Calculate saturation function
    integrator_arg <- ifelse(t-fit-1 > 0, sum(errs[1:(t-fit-1)]) - (t-fit-1)*alpha, 0)
    integrator <- saturation_fn_log(integrator_arg, t-fit, Csat, KI)
    
    # Train scorecaster if necessary
    if (scorecast & ((t-nfit-1) > nburnin)) {
      curr_scores <- scores[1:(t-nfit-1)]
      model <- forecast::thetaf(curr_scores, h = 1)
      scorecasts[t-nfit] <- as.numeric(model$mean)
    }
    
    # Update the next quantile
    if (t < T) {
      grad <- ifelse(errs[t-nfit], -(1-alpha), alpha)
      qts[t-nfit+1] <- qts[t-nfit] - lr_t * grad
      integrators[t-nfit+1] <- ifelse(integrate, integrator, 0)
      qs[t-nfit+1] <- qts[t-nfit+1] + integrators[t-nfit+1]
      if (scorecast) qs[t-nfit+1] <- qs[t-nfit+1] + scorecasts[t-nfit+1]
    }
  }
  
  return(qs)
}

saturation_fn_log <- function(x, t, Csat, KI) {
  if (KI == 0) {
    return(0)
  } else {
    tan_out <- mytan(x * log(t)/(Csat * t))
    out <- KI * tan_out
    return(out)
  }
}

mytan <- function(x){
  if (x >= pi/2) {
    return(Inf)
  } else if (x <= - pi/2) {
    return(-Inf)
  } else {
    return(tan(x))
  }
}
