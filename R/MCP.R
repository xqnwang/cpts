# Conformal forecasting using Multistep-ahead conformal prediction (MCP)
# This method can only deal with asymmetric conformity scores, i.e., forecast errors.
MCP <- function(object, alpha = 1 - 0.01 * object$level,
                ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = TRUE,
                lr = 0.1, Tg = length(object$x),
                Csat = 2 / pi * (ceiling(log(Tg) * 0.01) - 1 / log(Tg)),
                KI = abs(object$errors) |> max(na.rm = TRUE), ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  
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
  zeromatrix <- ts(matrix(0, nrow = NROW(errors), ncol = horizon),
                   start = start(errors),
                   frequency = frequency(errors))
  colnames(zeromatrix) <- paste0("h=", seq(horizon))
  
  lower <- upper <- rep(list(namatrix), length(alpha))
  names(lower) <- names(upper) <-
    paste0(level, "%")
  errt_lower <- errt_upper <-
    integrator_lower <- integrator_upper <- lower
  lrt <- scorecaster_lower <- scorecaster_upper <- namatrix
  
  out <- list(
    x = object$x,
    mean = object$mean,
    errors = object$errors,
    lower = lower,
    upper = upper,
    level = level
  )
  
  first_non_na <- min(which(!is.na(errors[, 1])))
  last_non_na <- max(which(!is.na(errors[, 1])))
  if (last_non_na < first_non_na + ncal - 1L)
    stop("errors in the input object is not long enough for calibration")
  indx <- seq(first_non_na + ncal - 1L, last_non_na, by = 1L)
  
  for (i in seq(length(alpha))) {
    lbl <- paste0(level[i], "%")
    qt_lower <- qt_upper <- zeromatrix
    
    for (t in indx) {
      errors_subset <- subset(
        errors,
        start = ifelse(!rolling, first_non_na, t - ncal + 1L),
        end = t)
      
      for (h in seq(horizon)) {
        # Calculate err_t
        errt_lower[[lbl]][t, h] <- (-errors[t, h] > qt_lower[t, h])
        errt_upper[[lbl]][t, h] <- (errors[t, h] > qt_upper[t, h])
        
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
        
        # Train scorecaster (same for different alphas)
        if (scorecast) {
          if (h == 1) {
            model <- forecast::meanf(errors_subset[, h], h = 1)
            scorecaster_upper[t, h] <- as.numeric(model$mean)
            scorecaster_lower[t, h] <- -scorecaster_upper[t, h]
            
          } else {
            model_MA <- forecast::Arima(errors_subset[, h], order = c(0, 0, h-1)) |>
              forecast::forecast(h = 1)
            model_LR <- lm(as.formula(paste0("V", h, " ~ .")), 
                           data = setNames(as.data.frame(errors_subset[, 1:h]), 
                                           paste0("V", 1:h))) |>
              forecast::forecast(newdata = setNames(data.frame(scorecaster_upper[t, 1:(h-1)] |> 
                                                                 matrix(nrow = 1)),
                                                    paste0("V", 1:(h-1))))
            scorecaster_upper[t, h] <- (as.numeric(model_MA$mean) + as.numeric(model_LR$mean))/2
            scorecaster_lower[t, h] <- -scorecaster_upper[t, h]
          }
        }
        
        # Learning rate (same for the upper and lower bounds)
        if (i == 1)
          lrt[t, h] <- ifelse(t == indx[1], lr,
                              lr*(max(errors_subset[, h]) - min(errors_subset[, h])))
        
        # Update the next quantile
        grad_lower <- ifelse(errt_lower[[lbl]][t, h], -(1-alpha[i]/2), alpha[i]/2)
        grad_upper <- ifelse(errt_upper[[lbl]][t, h], -(1-alpha[i]/2), alpha[i]/2)
        qt_lower[t+1, h] <- qt_lower[t, h] - lrt[t, h] * grad_lower +
          ifelse(integrate, integrator_lower[[lbl]][t, h], 0) +
          ifelse(scorecast, scorecaster_lower[t, h], 0)
        qt_upper[t+1, h] <- qt_upper[t, h] - lrt[t, h] * grad_upper +
          ifelse(integrate, integrator_upper[[lbl]][t, h], 0) +
          ifelse(scorecast, scorecaster_upper[t, h], 0)
        # PIs
        out$lower[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] - qt_lower[t+1, h]
        out$upper[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] + qt_upper[t+1, h]
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
  out$method <- paste("MCP")
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
