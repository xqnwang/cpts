# Conformal forecasting using adaptive conformal prediction method
ACP <- function(object, alpha = 1 - 0.01 * object$level, gamma = 0.005,
                symmetric = FALSE, ncal = 10, rollingwindow = FALSE,
                quantiletype = 1, ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  
  if (!quantiletype %in% 1:9)
    stop("quantiletype is invalid. It must be in 1:9.")
  
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  
  errors <- ts(as.matrix(object$errors),
               start = start(object$errors),
               frequency = frequency(object$errors))
  horizon <- NCOL(errors)
  
  namatrix <- ts(matrix(NA_real_, nrow = NROW(errors), ncol = horizon),
                 start = start(errors),
                 frequency = frequency(errors))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  lower <- upper <- rep(list(namatrix), length(alpha))
  names(lower) <- names(upper) <- paste0(level, "%")
  alphat <- alphat_lower <- alphat_upper <- rep(list(namatrix), length(alpha))
  names(alphat) <- names(alphat_lower) <- names(alphat_upper) <- paste0(level, "%")
  errt <- errt_lower <- errt_upper <- rep(list(namatrix), length(alpha))
  names(errt) <- names(errt_lower) <- names(errt_upper) <- paste0(level, "%")
  
  out <- list(
    x = object$x,
    mean = object$mean,
    errors = object$errors,
    lower = lower,
    upper = upper,
    level = level
  )
  
  for (h in seq(horizon)) {
    first_non_na <- min(which(!is.na(errors[, h])))
    last_non_na <- max(which(!is.na(errors[, h])))
    if (last_non_na < first_non_na + ncal - 1L)
      stop("errors in the input object is not long enough for calibration")
    indx <- seq(first_non_na + ncal - 1L, last_non_na, by = 1L)
    
    for (i in seq(length(alpha))) {
      lbl <- paste0(level[i], "%")
      if (symmetric) {
        alphat[[lbl]][indx[1]+1, h] <- alpha[i]
      } else {
        alphat_lower[[lbl]][indx[1]+1, h] <- alphat_upper[[lbl]][indx[1]+1, h] <-
          alpha[i]/2
      }
      
      for (t in indx) {
        errors_subset <- subset(
          errors[, h],
          start = ifelse(!rollingwindow, first_non_na, t - ncal + 1L),
          end = t)
        
        if (symmetric) {
          q_lo <- q_up <- ggdist::weighted_quantile(x = abs(c(errors_subset, Inf)),
                                                    probs = 1 - alphat[[lbl]][t+1, h],
                                                    type = quantiletype)
          # Compute errt
          if (alphat[[lbl]][t+1, h] >= 1) {
            errt[[lbl]][t+1, h] <- 1
          } else if (alphat[[lbl]][t+1, h] <= 0) {
            errt[[lbl]][t+1, h] <- 0
          } else {
            errt[[lbl]][t+1, h] <- as.numeric(abs(errors[t+1, h]) > q_lo)
          }
          # Update alphat
          alphat[[lbl]][t+2, h] <- alphat[[lbl]][t+1, h] +
            gamma*(alpha[i]- errt[[lbl]][t+1, h])
          
        } else {
          q_lo <- ggdist::weighted_quantile(x = - c(errors_subset, Inf),
                                            probs = 1 - alphat_lower[[lbl]][t+1, h],
                                            type = quantiletype)
          q_up <- ggdist::weighted_quantile(x = c(errors_subset, Inf),
                                            probs = 1 - alphat_upper[[lbl]][t+1, h],
                                            type = quantiletype)
          # Compute errt
          if (alphat_lower[[lbl]][t+1, h] >= 1) {
            errt_lower[[lbl]][t+1, h] <- 1
          } else if (alphat_lower[[lbl]][t+1, h] <= 0) {
            errt_lower[[lbl]][t+1, h] <- 0
          } else {
            errt_lower[[lbl]][t+1, h] <- as.numeric(-errors[t+1, h] > q_lo)
          }
          if (alphat_upper[[lbl]][t+1, h] >= 1) {
            errt_upper[[lbl]][t+1, h] <- 1
          } else if (alphat_lower[[lbl]][t+1, h] <= 0) {
            errt_upper[[lbl]][t+1, h] <- 0
          } else {
            errt_upper[[lbl]][t+1, h] <- as.numeric(errors[t+1, h] > q_up)
          }
          # Update alphat
          alphat_lower[[lbl]][t+2, h] <- alphat_lower[[lbl]][t+1, h] +
            gamma*(alpha[i]/2- errt_lower[[lbl]][t+1, h])
          alphat_upper[[lbl]][t+2, h] <- alphat_upper[[lbl]][t+1, h] +
            gamma*(alpha[i]/2- errt_upper[[lbl]][t+1, h])
        }
        out$lower[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] - q_lo
        out$upper[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] + q_up
      }
    }
  }
  if (h == 1) {
    out$lower <- lapply(out$lower, function(lo) lo[, 1L])
    out$upper <- lapply(upper$lower, function(up) up[, 1L])
    alphat <- lapply(alphat, function(alp) alp[, 1L])
    alphat_lower <- lapply(alphat_lower, function(alp_lo) alp_lo[, 1L])
    alphat_upper <- lapply(alphat_upper, function(alp_up) alp_up[, 1L])
  }
  out$method <- paste("ACP")
  if (symmetric) {
    out$model$alphat <- alphat
  } else {
    out$model$alphat_lower <- alphat_lower
    out$model$alphat_upper <- alphat_upper
  }
  return(structure(out, class = "CPforecast"))
}
