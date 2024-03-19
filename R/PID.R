# Conformal forecasting using conformal PID control
PID <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = !symmetric, scorecastfun = NULL,
                lr = 0.1, Tg = length(object$x),
                Csat = 2 / pi * (ceiling(log(Tg) * 0.01) - 1 / log(Tg)),
                KI = abs(object$errors) |> max(na.rm = TRUE), ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  if (scorecast && is.null(scorecastfun))
    stop("scorecastfun should not be NULL if scorecast is TRUE")
  
  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  pf <- ts(as.matrix(object$mean),
           start = start(object$mean),
           frequency = frequency(object$mean))
  errors <- ts(as.matrix(object$errors),
               start = start(object$errors),
               frequency = frequency(object$errors))
  horizon <- NCOL(errors)
  
  namatrix <- ts(matrix(NA_real_, nrow = NROW(errors), ncol = horizon), 
                 start = start(errors), 
                 frequency = frequency(errors))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  nalist <- rep(list(namatrix), length(alpha))
  names(nalist) <- paste0(level, "%")
  
  lower <- upper <- nalist
  lrmat <- namatrix
  if (integrate)
    integrator <- integrator_lower <- integrator_upper <- nalist
  if (scorecast)
    scorecaster <- scorecaster_lower <- scorecaster_upper <- namatrix
  
  out <- list(
    x = object$x
  )
  
  for (h in seq(horizon)) {
    first_non_na <- (!is.na(errors[, h])) |> which() |> min()
    last_non_na <- (!is.na(errors[, h])) |> which() |> max()
    if (last_non_na < first_non_na + ncal - 1L)
      stop("errors in the input object is not long enough for calibration")
    indx <- seq(first_non_na + ncal - 1L, last_non_na + h - 1L, by = 1L)
    
    errt_h <- errt_lower_h <- errt_upper_h <-
      integ_h <- integ_lower_h <- integ_upper_h <-
      matrix(NA_real_, nrow = NROW(errors), ncol = length(alpha))
    qts_h <- qts_lower_h <- qts_upper_h <- matrix(0, nrow = NROW(errors), ncol = length(alpha))
    qs_h <- qs_lower_h <- qs_upper_h <- matrix(0, nrow = NROW(errors), ncol = length(alpha))
    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, first_non_na, t - ncal + 1L),
        end = ifelse(t <= last_non_na, t, last_non_na))
      
      if (symmetric) {
        if (t <= last_non_na) {
          # Calculate errt
          errt_h[t, ] <- (abs(errors[t, h]) > qs_h[t, ])
          
          # Learning rate (same for the upper and lower bounds)
          lrmat[t, h] <- ifelse(
            length(errors_subset) <= 1,
            lr,
            lr*(max(abs(errors_subset)) - min(abs(errors_subset))))
          
          # Update quantile tracking
          qts_h[t+1, ] <- qts_h[t, ] + lrmat[t, h] * (errt_h[t, ] - alpha)
        } else {
          qts_h[t+1, ] <- qts_h[t, ]
        }
        
        # Update integrator
        if (integrate) {
          es <- errt_h[indx[1]:(ifelse(t <= last_non_na, t, last_non_na)), ] |>
            matrix(ncol = length(alpha))
          integrator_arg <- apply(es, 2, sum) - NROW(es)*alpha
          integ_h[t+1, ] <- sapply(
            1:length(alpha),
            function(i) ifelse(
              NROW(es) == 1,
              0,
              saturation_fn_log(integrator_arg[i], NROW(es), Csat, KI)))
        }
        
        # Update scorecaster
        if (scorecast) {
          sc <- try(suppressWarnings(
            scorecastfun(abs(errors_subset), h = 1, ...)
          ), silent = TRUE)
          if (!is.element("try-error", class(sc))) {
            scorecaster[t+1, h] <- as.numeric(sc$mean, 1)
          }
        }
        
        # Update the next quantile
        qs_h[t+1, ] <- qts_h[t+1, ] +
          ifelse(rep(integrate, length(alpha)), integ_h[t+1, ], rep(0, length(alpha))) +
          rep(ifelse(scorecast, scorecaster[t+1, h], 0), length(alpha))
        qs_lower_h[t+1, ] <- qs_upper_h[t+1, ] <- qs_h[t+1, ]
      } else {
        if (t <= last_non_na) {
          # Calculate errt
          errt_lower_h[t, ] <- (-errors[t, h]) > qs_lower_h[t, ]
          errt_upper_h[t, ] <- errors[t, h] > qs_upper_h[t, ]
          
          # Learning rate (same for the upper and lower bounds)
          lrmat[t, h] <- ifelse(
            length(errors_subset) <= 1,
            lr,
            lr*(max(errors_subset) - min(errors_subset)))
          
          # Update quantile tracking
          qts_lower_h[t+1, ] <- qts_lower_h[t, ] + lrmat[t, h] * (errt_lower_h[t, ] - alpha/2)
          qts_upper_h[t+1, ] <- qts_upper_h[t, ] + lrmat[t, h] * (errt_upper_h[t, ] - alpha/2)
        } else {
          qts_lower_h[t+1, ] <- qts_lower_h[t, ]
          qts_upper_h[t+1, ] <- qts_upper_h[t, ]
        }
        
        # Update integrator
        if (integrate) {
          el <- errt_lower_h[indx[1]:(ifelse(t <= last_non_na, t, last_non_na)), ] |>
            matrix(ncol = length(alpha))
          integrator_lower_arg <- apply(el, 2, sum) - NROW(el)*alpha/2
          integ_lower_h[t+1, ] <- sapply(
            1:length(alpha),
            function(i) saturation_fn_log(integrator_lower_arg[i], NROW(el), Csat, KI))
          
          eu <- errt_upper_h[indx[1]:(ifelse(t <= last_non_na, t, last_non_na)), ] |>
            matrix(ncol = length(alpha))
          integrator_upper_arg <- apply(eu, 2, sum) - NROW(eu)*alpha/2
          integ_upper_h[t+1, ] <- sapply(
            1:length(alpha),
            function(i) saturation_fn_log(integrator_upper_arg[i], NROW(eu), Csat, KI))
        }
        
        # Update scorecaster
        if (scorecast) {
          sc <- try(suppressWarnings(
            scorecastfun(errors_subset, h = 1, ...)
          ), silent = TRUE)
          if (!is.element("try-error", class(sc))) {
            scorecaster_lower[t+1, h] <- - as.numeric(sc$mean)
            scorecaster_upper[t+1, h] <- as.numeric(sc$mean)
          }
        }
        
        # Update the next quantile
        qs_lower_h[t+1, ] <- qts_lower_h[t+1, ] +
          ifelse(rep(integrate, length(alpha)), integ_lower_h[t+1, ], rep(0, length(alpha))) +
          rep(ifelse(scorecast, scorecaster_lower[t+1, h], 0), length(alpha))
        qs_upper_h[t+1, ] <- qts_upper_h[t+1, ] +
          ifelse(rep(integrate, length(alpha)), integ_upper_h[t+1, ], rep(0, length(alpha))) +
          rep(ifelse(scorecast, scorecaster_upper[t+1, h], 0), length(alpha))
      }
      
      # PIs
      for (i in seq(length(alpha))) {
        lbl <- paste0(level[i], "%")
        lower[[lbl]][t+1, h] <- pf[t+1, h] - qs_lower_h[t+1, i]
        upper[[lbl]][t+1, h] <- pf[t+1, h] + qs_upper_h[t+1, i]
      }
    }
    if (integrate) {
      for (i in seq(length(alpha))) {
        integrator[[i]][, h] <- integ_h[, i]
        integrator_lower[[i]][, h] <- integ_lower_h[, i]
        integrator_upper[[i]][, h] <- integ_upper_h[, i]
      }
    }
  }
  if (h == 1) {
    lower <- lapply(out$lower, function(lo) lo[, 1L])
    upper <- lapply(out$upper, function(up) up[, 1L])
    lrmat <- lrmat[, 1L]
    if (integrate) {
      integrator <- lapply(integrator, function(integrator_t) integrator_t[, 1L])
      integrator_lower <- lapply(integrator_lower, function(int_low) int_low[, 1L])
      integrator_upper <- lapply(integrator_upper, function(int_up) int_up[, 1L])
    }
    if (scorecast) {
      scorecaster <- scorecaster[, 1L]
      scorecaster_lower <- scorecaster_lower[, 1L]
      scorecaster_upper <- scorecaster_upper[, 1L]
    }
  }
  out$method <- paste("PID")
  out$mean <- object$mean
  out$errors <- object$errors
  out$lower <- lower
  out$upper <- upper
  out$level <- level
  out$model$method <- out$method
  out$model$call <- match.call()
  out$model$alpha <- alpha
  out$model$symmetric <- symmetric
  out$model$ncal <- ncal
  out$model$rolling <- rolling
  out$model$integrate <- integrate
  out$model$scorecast <- scorecast
  out$model$lr <- lr
  out$model$Tg <- Tg
  out$model$Csat <- Csat
  out$model$KI <- KI
  if (symmetric) {
    if (integrate)
      out$model$integrator <- integrator
    if (scorecast)
      out$model$scorecaster <- scorecaster
  } else {
    if (integrate)
      out$model$integrator <- list(lower = integrator_lower, upper = integrator_upper)
    if (scorecast)
      out$model$scorecaster <- list(lower = scorecaster_lower, upper = scorecaster_upper)
  }
  
  return(structure(out, class = "CPforecast"))
}

saturation_fn_log <- function(x, t, Csat, KI) {
  if (KI == 0) {
    return(0)
  } else {
    tan_out <- mytan(x * log(t)/(Csat * (t)))
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