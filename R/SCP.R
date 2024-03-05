# Conformal forecasting using classical split conformal prediction method
SCP <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rollingwindow = FALSE,
                quantiletype = 1, weightfunction = NULL, kess = FALSE,...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  
  if (is.null(weightfunction)) {
    weightfunction <- function(n) rep(1, n + 1L)
  }
  
  # Kish's effective sample size for sample quantile computation
  if (kess) {
    kess <- function(w) sum(w)^2 / sum(w^2)
  } else {
    kess <- NULL
  }
  
  if (!quantiletype %in% 1:9)
    stop("quantiletype is invalid. It must be in 1:9.")
  
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
    
    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rollingwindow, first_non_na, t - ncal + 1L),
        end = t)
      
      weight_subset <- weightfunction(length(errors_subset) + 1L, ...)
      
      for (i in seq(length(alpha))) {
        if (symmetric) {
          q_lo <- q_up <- ggdist::weighted_quantile(x = abs(c(errors_subset, Inf)),
                                                    probs = 1 - alpha[i],
                                                    weights = weight_subset,
                                                    n = kess,
                                                    type = quantiletype)
        } else {
          q_lo <- ggdist::weighted_quantile(x = - c(errors_subset, Inf),
                                            probs = 1 - alpha[i]/2,
                                            weights = weight_subset,
                                            n = kess,
                                            type = quantiletype)
          q_up <- ggdist::weighted_quantile(x = c(errors_subset, Inf),
                                            probs = 1 - alpha[i]/2,
                                            weights = weight_subset,
                                            n = kess,
                                            type = quantiletype)
        }
        out$lower[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] - q_lo
        out$upper[[paste0(level[i], "%")]][t+1, h] <- out$mean[t+1, h] + q_up
      }
    }
  }
  if (h == 1) {
    out$lower <- lapply(out$lower, function(lo) lo[, 1L])
    out$upper <- lapply(upper$lower, function(up) up[, 1L])
  }
  out$method <- paste("SCP")
  return(structure(out, class = "CPforecast"))
}
