PID.update <- function(series, nfit, nburnin = 100, alpha = 0.1,
                       integrate = TRUE, scorecast = TRUE, ncast = NULL,
                       lr = 0.1, Csat = NULL, KI = NULL) {
  # Set up data
  y <- series
  T <- length(y)
  
  # Check the series length
  if (T < (nfit + nburnin + 1))
    stop("series length must be larger than nfit+nburnin")
  
  # Initialize data storage variables
  predSeq <- loSeq <- upSeq <- rep(0, T - nfit - nburnin)
  scores <- errs <- qs <- qts <- integrators <- scorecasts <- rep(0, T - nfit)
  
  for (t in (nfit + 1):T) {
    # Learning rate
    t_lr <- t - nfit
    t_lr_min <- max(t_lr - nburnin, 1)
    lr_t <- ifelse(t_lr > 1, lr * (max(scores[t_lr_min:(t_lr - 1)]) - min(scores[t_lr_min:(t_lr - 1)])), lr)
    
    # Fit AR model and get point forecasts
    out <- forecast::Arima(y[t_lr:(t-1)], # use rolling window (length = nfit) for AR model training
                           order = c(2,0,0),
                           include.mean = TRUE, method = "CSS")
    fc <- forecast::forecast(out, h = 1)
    pred <- as.vector(fc$mean)
    
    # Observe y_t, calculate conformity score and err
    scores[t_lr] <- abs(y[t] - pred)
    errs[t_lr] <- !(qs[t_lr] >= scores[t_lr])
    
    # Calculate saturation function
    integrator_arg <- ifelse(t_lr - 1 > 0, sum(errs[1:(t_lr - 1)]) - (t_lr - 1) * alpha, 0)
    integrator <- saturation_fn_log(integrator_arg, t_lr, Csat, KI)
    
    # Train scorecaster if necessary
    if (scorecast && (nburnin <= t_lr) && (t_lr < (T - nfit))) {
      if (is.null(ncast) | (ncast > nburnin)) {
        curr_scores <- scores[1:t_lr] # use expanding window for scorecaster
      } else {
        curr_scores <- scores[(t_lr - ncast + 1):t_lr] # use rolling window for scorecaster, normally set it to nburnin
      }
      model <- forecast::thetaf(curr_scores, h = 1)
      scorecasts[t_lr + 1] <- as.numeric(model$mean)
    }
    
    # Update the next quantile
    if (t_lr <= (T - nfit - 1)) {
      grad <- ifelse(errs[t_lr], -(1-alpha), alpha)
      qts[t_lr + 1] <- qts[t_lr] - lr_t * grad # set q_1 = 0
      integrators[t_lr + 1] <- ifelse(integrate, integrator, 0)
      qs[t_lr + 1] <- qts[t_lr + 1] + integrators[t_lr + 1]
      if (scorecast) qs[t_lr + 1] <- qs[t_lr + 1] + scorecasts[t_lr + 1]
    }
    
    # Get forecasts and PIs
    if (t_lr > nburnin) {
      lo <- pred - qs[t_lr]
      up <- pred + qs[t_lr]
      
      predSeq[t_lr - nburnin] <- pred
      loSeq[t_lr - nburnin] <- min(lo, up)
      upSeq[t_lr - nburnin] <- max(lo, up)
    }
  }
  return(list(pred = predSeq, lo = loSeq, up = upSeq))
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
