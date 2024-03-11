# Conformal forecasting using conformal PID control
PID <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = !symmetric, scorecastfun = NULL,
                lr = 0.1, Csat = NULL, KI = NULL, ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  
  if (scorecast) {
    if (is.null(scorecastfun))
      stop("scorecastfun should not be NULL if scorecast is TRUE")
  }
  
  errors <- ts(as.matrix(object$errors),
               start = start(object$errors),
               frequency = frequency(object$errors))
  horizon <- NCOL(errors)
  
  namatrix <- ts(matrix(NA_real_, nrow = NROW(errors), ncol = horizon), 
                 start = start(errors), 
                 frequency = frequency(errors))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  lrt <- namatrix
  lower <- upper <- rep(list(namatrix), length(alpha))
  names(lower) <- names(upper) <- paste0(level, "%")
  if (symmetric) {
    errt <- integrator <- lower
    scorecaster <- namatrix
  } else {
    errt_lower <- errt_upper <-
      integrator_lower <- integrator_upper <- lower
    scorecaster_lower <- scorecaster_upper <- namatrix
  }
  
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
        qs <- 0
      } else {
        qs_lower <- qs_upper <- 0
      }
      
      for (t in indx) {
        errors_subset <- subset(
          errors[, h],
          start = ifelse(!rolling, first_non_na, t - ncal + 1L),
          end = t)
        
        if (symmetric) {
          # Calculate err_t
          errt[[lbl]][t, h] <- as.numeric(abs(errors[t, h]) > qs)
          # Calculate saturation function using err_i, i = indx[1],...,t-1
          if (integrate) {
            integrator_arg <- ifelse(t > indx[1],
                                     sum(errt[[lbl]][indx[1]:(t-1), h]) - (t - indx[1]) * alpha[i],
                                     0)
            integrator[[lbl]][t, h] <- saturation_fn_log(integrator_arg, t-indx[1]+1, Csat, KI)
          }
          # Train scorecaster
          if (scorecast && i == 1) {
            sc <- try(suppressWarnings(
              scorecastfun(abs(errors_subset), h = 1, ...)
            ), silent = TRUE)
            if (!is.element("try-error", class(sc))) {
              scorecaster[t, h] <- as.numeric(sc$mean)
            }
          }
          # Learning rate
          if (i == 1)
            lrt[t, h] <- ifelse(t == indx[1], lr,
                                lr*(max(abs(errors_subset)) - min(abs(errors_subset))))
          # Update the next quantile
          grad <- ifelse(errt[[lbl]][t, h], -(1-alpha[i]), alpha[i])
          qs <- qs - lrt[t, h] * grad +
            ifelse(integrate, integrator[[lbl]][t, h], 0) +
            ifelse(scorecast, scorecaster[t, h], 0)
          # PIs
          out$lower[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] - qs
          out$upper[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] + qs
          
        } else {
          # Calculate err_t
          errt_lower[[lbl]][t, h] <- as.numeric(-errors[t, h] > qs_lower)
          errt_upper[[lbl]][t, h] <- as.numeric(errors[t, h] > qs_upper)
          # Calculate saturation function using err_i, i = indx[1],...,t-1
          if (integrate) {
            integrator_arg_lower <- ifelse(t > indx[1],
                                           sum(errt_lower[[lbl]][indx[1]:(t-1), h]) -
                                             (t - indx[1]) * alpha[i]/2,
                                           0)
            integrator_arg_upper <- ifelse(t > indx[1],
                                           sum(errt_upper[[lbl]][indx[1]:(t-1), h]) -
                                             (t - indx[1]) * alpha[i]/2,
                                           0)
            integrator_lower[[lbl]][t, h] <- saturation_fn_log(integrator_arg_lower, t-indx[1]+1, Csat, KI)
            integrator_upper[[lbl]][t, h] <- saturation_fn_log(integrator_arg_upper, t-indx[1]+1, Csat, KI)
          }
          # Train scorecaster
          if (scorecast && i == 1) {
            sc_upper <- try(suppressWarnings(
              scorecastfun(errors_subset, h = 1, ...)
            ), silent = TRUE)
            if (!is.element("try-error", class(sc_upper))) {
              scorecaster_lower[t, h] <- - as.numeric(sc_upper$mean)
              scorecaster_upper[t, h] <- as.numeric(sc_upper$mean)
            }
          }
          # Learning rate (same for the upper and lower bounds)
          if (i == 1)
            lrt[t, h] <- ifelse(t == indx[1], lr,
                                lr*(max(errors_subset) - min(errors_subset)))
          # Update the next quantile
          grad_lower <- ifelse(errt_lower[[lbl]][t, h], -(1-alpha[i]/2), alpha[i]/2)
          grad_upper <- ifelse(errt_upper[[lbl]][t, h], -(1-alpha[i]/2), alpha[i]/2)
          qs_lower <- qs_lower - lrt[t, h] * grad_lower +
            ifelse(integrate, integrator_lower[[lbl]][t, h], 0) +
            ifelse(scorecast, scorecaster_lower[t, h], 0)
          qs_upper <- qs_upper - lrt[t, h] * grad_upper +
            ifelse(integrate, integrator_upper[[lbl]][t, h], 0) +
            ifelse(scorecast, scorecaster_upper[t, h], 0)
          # PIs
          out$lower[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] - qs_lower
          out$upper[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] + qs_upper
        }
      }
    }
  }
  if (h == 1) {
    out$lower <- lapply(out$lower, function(lo) lo[, 1L])
    out$upper <- lapply(upper$lower, function(up) up[, 1L])
    lrt <- lrt[, 1L]
    if (symmetric) {
      integrator <- lapply(integrator, function(integrator_t) integrator_t[, 1L])
      scorecaster <- scorecaster[, 1L]
    } else {
      integrator_lower <- lapply(integrator_lower, function(integrator_lower_t) integrator_lower_t[, 1L])
      integrator_upper <- lapply(integrator_upper, function(integrator_upper_t) integrator_upper_t[, 1L])
      scorecaster_lower <- scorecaster_lower[, 1L]
      scorecaster_upper <- scorecaster_upper[, 1L]
    }
  }
  out$method <- paste("PID")
  out$model$lr <- lrt
  if (symmetric) {
    out$model$integrator <- list(integrator = integrator)
    out$model$scorecaster <- list(scorecaster = scorecaster)
  } else {
    out$model$integrator <- list(integrator_lower = integrator_lower, integrator_upper = integrator_upper)
    out$model$scorecaster <- list(scorecaster_lower = scorecaster_lower, scorecaster_upper = scorecaster_upper)
  }
  return(structure(out, class = "CPforecast"))
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