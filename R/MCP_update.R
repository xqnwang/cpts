MCP.update <- function(series, h = 1, nfit, nburnin = 100, alpha = 0.1,
                       integrate = TRUE, scorecast = TRUE, ncast = NULL,
                       lr = 0.1, Csat, KI) {
  # Set up data
  y <- series
  T <- length(y)
  
  # Check the series length
  if (T < (nfit + nburnin + 1))
    stop("series length must be larger than nfit+nburnin")
  
  # Initialize data storage variables
  predSeq <- loSeq <- upSeq <- matrix(0, nrow = T - nfit - nburnin - h + 1, ncol = h)
  
  for (upper in c(FALSE, TRUE)) {
    # Initialize data storage variables
    scores <- errs <- qs <- qts <- integrators <- scorecasts <- matrix(0, nrow = T - nfit - h + 1, ncol = h)
    
    for (t in (nfit + 1):(T - h + 1)) {
      t_lr <- t - nfit
      t_lr_min <- max(t_lr - nburnin, 1)
      
      # Fit AR model and get point forecasts
      out <- forecast::Arima(y[t_lr:(t-1)], # use rolling window (length = nfit) for AR model training
                             order = c(2,0,0),
                             include.mean = TRUE, method = "CSS")
      fc <- forecast::forecast(out, h = h)
      pred <- as.vector(fc$mean)
      
      # Observe y_t, calculate conformity score and err
      if (upper) {
        scores[t_lr, ] <- y[t:(t+h-1)] - pred
      } else {
        scores[t_lr, ] <- pred - y[t:(t+h-1)]
      }
      errs[t_lr, ] <- !(qs[t_lr, ] >= scores[t_lr, ])
      
      for (j in 1:h) {
        # Calculate saturation function
        integrator_arg <- ifelse(t_lr - 1 > 0, sum(errs[1:(t_lr - 1), j]) - (t_lr - 1) * alpha/2, 0)
        integrator <- saturation_fn_log(integrator_arg, t_lr, Csat, KI)
        
        # Train scorecaster if necessary
        if (scorecast && (nburnin <= t_lr) && (t_lr < (T - nfit - h + 1))) {
          if (is.null(ncast) | (ncast > nburnin)) {
            curr_scores <- scores[1:t_lr, ] # use expanding window for scorecaster
          } else {
            curr_scores <- scores[(t_lr - ncast + 1):t_lr, ] # use rolling window for scorecaster, normally set it to nburnin
          }
          
          if (j == 1){
            model <- forecast::meanf(curr_scores[, j], # use rolling window (length = nfit) for AR model training
                                      h = 1)
            
            scorecasts[t_lr + 1, j] <- as.numeric(model$mean)
          } else {
            model_MA <- forecast::Arima(curr_scores[, j], order = c(0, 0, j-1),
                                        include.mean = TRUE, method = "CSS") |>
              forecast::forecast(h = 1)
            model_LR <- lm(as.formula(paste0("V", j, " ~ .")), 
                           data = setNames(as.data.frame(curr_scores[, 1:j]), 
                                           paste0("V", 1:j))) |>
              forecast::forecast(newdata = setNames(data.frame(scorecasts[t_lr + 1, 1:(j-1)] |> matrix(nrow=1)),
                                                    paste0("V", 1:(j-1))))
            
            scorecasts[t_lr + 1, j] <- (as.numeric(model_MA$mean) + as.numeric(model_LR$mean))/2
          }
        }
        
        # Learning rate
        lr_t <- ifelse(t_lr > 1, lr * (max(scores[t_lr_min:(t_lr - 1), j]) - min(scores[t_lr_min:(t_lr - 1), j])), lr)
        
        # Update the next quantile
        if (t_lr <= (T - nfit - h)) {
          grad <- ifelse(errs[t_lr, j], -(1-alpha/2), alpha/2)
          qts[t_lr + 1, j] <- qts[t_lr, j] - lr_t * grad # set q_1 = 0
          integrators[t_lr + 1, j] <- ifelse(integrate, integrator, 0)
          qs[t_lr + 1, j] <- qts[t_lr + 1, j] + integrators[t_lr + 1, j]
          if (scorecast) qs[t_lr + 1, j] <- qs[t_lr + 1, j] + scorecasts[t_lr + 1, j]
        }
      }
      
      # Get forecasts and PIs
      if (t_lr > nburnin) {
        predSeq[t_lr - nburnin, ] <- pred
        if (upper) {
          upSeq[t_lr - nburnin, ] <- pred + qs[t_lr, ]
        } else {
          loSeq[t_lr - nburnin, ] <- pred - qs[t_lr, ]
        }
      }
    }
  }
  return(list(pred = predSeq, lo = loSeq, up = upSeq))
}
