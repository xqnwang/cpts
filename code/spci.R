# Multi-step SPCI using rolling out-of-sample signed forecast errors.
# Requires reticulate and Python package sklearn-quantile.
# w is the residual lag length; ... is passed to RandomForestQuantileRegressor().
spci <- function(object, alpha = 1 - 0.01 * object$level,
                 ncal = 10, rolling = FALSE, w = 5,
                 na.rm = TRUE, num.trees = 10, max.depth = 2,
                 mtry = w, min.node.size = 1, ...) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package `reticulate` is required")
  qrf <- reticulate::import("sklearn_quantile")
  if (!length(alpha) || anyNA(alpha) || any(alpha >= 1 | alpha <= 0))
    stop("`alpha` should be in (0, 1)")
  if (length(ncal) != 1 || is.na(ncal) || !is.finite(ncal) ||
      ncal < 10 || ncal != floor(ncal))
    stop("`ncal` should be an integer of at least 10")
  if (length(w) != 1 || is.na(w) || !is.finite(w) ||
      w < 1 || w != floor(w))
    stop("`w` should be a positive integer")

  alpha <- sort(alpha, decreasing = TRUE)
  level <- 100 * (1 - alpha)
  pf <- ts(as.matrix(object$MEAN), start = start(object$MEAN),
           frequency = frequency(object$MEAN))
  errors <- ts(as.matrix(object$ERROR), start = start(object$ERROR),
               frequency = frequency(object$ERROR))
  horizon <- ncol(pf)
  n <- nrow(pf)
  if (ncal > nrow(errors))
    stop("`ncal` is larger than the number of rows in object$ERROR")
  if (ncal < w + horizon)
    stop("`ncal` should be at least `w` plus the forecast horizon")

  namatrix <- ts(matrix(NA_real_, nrow = n, ncol = horizon),
                 start = start(pf), frequency = frequency(pf))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  lower <- upper <- setNames(rep(list(namatrix), length(alpha)),
                             paste0(level, "%"))
  out <- c(list(x = object$x, series = object$series),
           if ("xreg" %in% names(object)) list(xreg = object$xreg))
  probs <- c(alpha / 2, 1 - alpha / 2)
  quantiles <- sort(unique(probs))

  for (h in seq(horizon)) {
    first <- ncal + h - 1L
    last <- nrow(errors) - !object$forward
    indx <- if (first <= last) seq.int(first, last) else integer(0)
    for (t in indx) {
      t_start <- if (rolling) t - ncal + 1L else h
      target <- seq.int(t_start + w + h - 1L, t)

      # At origin u-h, only errors up to u-h can predict error u.
      x <- vapply(seq_len(w), function(j) {
        errors[target - h - j + 1L, h]
      }, numeric(length(target)))
      x <- matrix(x, ncol = w)
      colnames(x) <- paste0("lag", seq_len(w))
      y <- as.numeric(errors[target, h])
      newx <- matrix(errors[t - seq_len(w) + 1L, h], nrow = 1)
      colnames(newx) <- colnames(x)
      keep <- complete.cases(x, y)
      if (!na.rm && (anyNA(newx) || !all(keep)))
        stop("Missing calibration errors; use `na.rm = TRUE`")
      if (anyNA(newx) || is.na(pf[t + h, h]) || !any(keep))
        next

      model <- qrf$RandomForestQuantileRegressor(
        n_estimators = as.integer(num.trees),
        max_depth = if (is.null(max.depth) || max.depth == 0) NULL else as.integer(max.depth),
        max_features = as.integer(mtry),
        min_samples_split = as.integer(min.node.size + 1L),
        min_samples_leaf = 1L, criterion = "squared_error",
        bootstrap = TRUE, n_jobs = 1L, q = quantiles, ...)
      model$fit(x[keep, , drop = FALSE], array(y[keep]))
      q <- as.numeric(model$predict(newx))[match(probs, quantiles)]
      for (i in seq_along(alpha)) {
        lower[[i]][t + h, h] <- pf[t + h, h] + q[i]
        upper[[i]][t + h, h] <- pf[t + h, h] + q[length(alpha) + i]
      }
    }
  }

  out$method <- "spci"
  out$cp_times <- length(indx)
  out$MEAN <- object$MEAN
  out$ERROR <- object$ERROR
  out$LOWER <- lower
  out$UPPER <- upper
  out$level <- level
  out$forward <- object$forward
  out$call <- match.call()
  if ("mean" %in% names(object)) {
    extract_final <- getFromNamespace("extract_final", "conformalForecast")
    out$mean <- object$mean
    out$lower <- extract_final(lower, nrow = n, ncol = horizon,
                               bench = out$mean)
    out$upper <- extract_final(upper, nrow = n, ncol = horizon,
                               bench = out$mean)
  }
  out$model$method <- out$method
  out$model$call <- match.call()
  out$model$alpha <- alpha
  out$model$ncal <- ncal
  out$model$rolling <- rolling
  out$model$w <- w
  out$model$num.trees <- num.trees
  out$model$max.depth <- max.depth
  out$model$mtry <- mtry
  out$model$min.node.size <- min.node.size
  out$model$cvforecast$call <- object$call
  return(structure(out, class = c("spci", "cpforecast", "forecast")))
}
