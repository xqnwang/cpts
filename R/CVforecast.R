CVforecast <- function(y, forecastfunction, h = 1, level = c(80, 95), PI = TRUE, forward = TRUE, window = NULL, xreg = NULL, initial = 0, ...) {
  y <- as.ts(y)
  n <- length(y)
  # Order levels
  level <- sort(level)
  if (h <= 0) 
    stop("Forecast horizon out of bounds")
  if (initial >= n)
    stop("initial period too long")
  if (forward) {
    N <- n + h
  } else {
    N <- n
  }
  if (!is.null(xreg)) {
    if (NROW(xreg) != N)
      stop("xreg must be of the same size as y (if forward = FALSE) or 
           the same size as y plus h (if forward = TRUE)")
    # Make xreg a ts object to allow easy subsetting later
    xreg <- ts(as.matrix(xreg), start = start(y), frequency = frequency(y))
  }
  if (is.null(window)) {
    indx <- seq(1 + initial, n, by = 1L)
  } else {
    indx <- seq(window + initial, n, by = 1L)
  }
  
  pf <- err <- ts(matrix(NA_real_, nrow = N, ncol = h), 
                  start = start(y), 
                  frequency = frequency(y))
  colnames(pf) <- colnames(err) <- paste0("h=", seq(h))
  if (PI) {
    lower <- rep(list(pf), length(level))
    upper <- rep(list(pf), length(level))
    names(lower) <- names(upper) <- paste0(level, "%")
  }
  out <- list(
    x = y
  )
  
  for (i in indx) {
    y_subset <- subset(
      y,
      start = ifelse(is.null(window), 1L,
                     ifelse(i - window >= 0L, i - window + 1L, stop("small window"))),
      end = i)
    
    if (is.null(xreg)) {
      fc <- try(suppressWarnings(
        forecastfunction(y_subset, h = h, level = level, PI = TRUE, ...)
        ), silent = TRUE)
    } else {
      xreg_subset <- subset(
        xreg,
        start = ifelse(is.null(window), 1L,
                       ifelse(i - window >= 0L, i - window + 1L, stop("small window"))),
        end = i)
      xreg_future <- subset(
        xreg,
        start = i + 1L,
        end = i + h)
      fc <- try(suppressWarnings(
        forecastfunction(y_subset, h = h, level = level, PI = TRUE, xreg = xreg_subset, newxreg = xreg_future, ...)
        ), silent = TRUE)
    }
    
    if (!is.element("try-error", class(fc))) {
      pf[i,] <- fc$mean[seq(h)]
      err[i,] <- y[i + seq(h)] - fc$mean[seq(h)]
      if (PI) {
        for (l in level) {
          lower[[paste0(l, "%")]][i,] <- fc$lower[seq(h), paste0(l, "%")]
          upper[[paste0(l, "%")]][i,] <- fc$upper[seq(h), paste0(l, "%")]
        }
      }
    }
  }
  
  out$mean <- lagmatrix(pf, seq(h))
  out$errors <- lagmatrix(err, seq(h))
  if (!PI) {
    out$lower <- out$upper <- out$level <- NULL
  } else {
    out$lower <- lapply(lower, function(low) lagmatrix(low, seq(h)))
    out$upper <- lapply(upper, function(up) lagmatrix(up, seq(h)))
    out$level <- level
  }
  if (!is.null(xreg)){
    out$xreg <- xreg
  }
  # The final forecasting model output in the for loop
  out$method <- fc$method
  out$model <- fc$model
  out$fitted <- fc$fitted
  out$series <- fc$series
  out$residuals <- fc$residuals
  
  if (h == 1) {
    out$mean <- out$mean[, 1L]
    out$errors <- out$errors[, 1L]
    out$lower <- out$lower[, 1L]
    out$upper <- out$upper[, 1L]
  }
  
  return(structure(out, class = "CVforecast"))
}

lagmatrix <- function(x, k) {
  # ensure 'x' is a matrix
  if (!is.matrix(x))
    stop("ensure x is a matrix")
  n <- NROW(x)
  if (NCOL(x) != length(k))
    stop("k must have the same number of columns as x")
  y <- x
  if (all(k == rep(0, length(k)))) {
    return(y)
  } else {
    for (i in seq(NCOL(x))) {
      if (k[i] == 0) {
        y[, i] <- x[, i]
      } else {
        y[, i] <- c(rep(NA, k[i]), x[-(seq(n - k[i] + 1, n, by = 1L)), i])
      }
    }
  }
  
  return(y)
}
