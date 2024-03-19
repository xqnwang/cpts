#' Adaptive conformal prediction methods for time series forecasting
#' 
#' \code{ACP} computes prediction intervals and other information obtained by
#' applying the adaptive conformal prediction method.
#' 
#' @param object An object of class "\code{CVforecast}". It must have an argument
#' \code{x} for original univariate time series, an argument \code{mean} for
#' point forecasts and \code{errors} for forecast errors. See the results of a call
#' to \code{\link{CVforecast}}.
#' @param alpha A numeric vector of target levels \eqn{1-\alpha}.
#' @param gamma The step size parameter \eqn{\gamma>0} for the \eqn{\alpha} update.
#' @param symmetric If \code{TRUE}, symmetric conformity scores (i.e. \eqn{|e_{t+h|t}|})
#' are used. If \code{FALSE}, asymmetric conformity scores are used, and upper
#' limits and lower limits are generated separately.
#' @param ncal Length of the calibration set. If \code{rolling=FALSE}, it denotes
#' the initial period of calibration sets. If \code{rolling=TRUE}, it indicates
#' the period of each rolling calibration set.
#' \code{rolling=TRUE}.
#' @param rolling If \code{TRUE}, a rolling window strategy will be used to
#' generate the calibration set. Otherwise, expanding window will be used.
#' @param quantiletype An integer between 1 and 9 determining the type of
#' quantile estimator to be used. Types 1 to 3 are for discontinuous quantiles,
#' types 4 to 9 are for continuous quantiles. See the
#' \code{\link[ggdist]{weighted_quantile}} function in the ggdist package.
#' @param na.rm If \code{TRUE}, \code{NA} values are removed when calculating
#' sample quantile.
#' @param ... Other arguments are passed to the
#' \code{\link[ggdist]{weighted_quantile}} function in the ggdist package for
#' sample quantile computation.
#' @return An object of class "\code{CPforecast}".
#' 
#' An object of class "\code{CPforecast}" is a list containing the following elements:
#' \item{x}{The original time series.}
#' \item{method}{The name of the conformal prediction method as a character string.}
#' \item{mean}{Point forecasts as a time series for \eqn{h=1}. For \eqn{h>1},
#' they are returned as a time series matrix with the \eqn{h}th column
#' containing point forecasts for forecast horizon \eqn{h}. The time index
#' corresponds to the period for which the forecast is generated.}
#' \item{errors}{Forecast errors given by
#' \eqn{e_{t+h} = y_{t+h}-\hat{y}_{t+h|t}}{e[t+h] = y[t+h]-f[t+h]}.}
#' \item{lower}{A list containing lower limits for prediction intervals for
#' each \code{level}.}
#' \item{upper}{A list containing upper limits for prediction intervals for
#' each \code{level}.}
#' \item{level}{The confidence values associated with the prediction intervals.}
#' \item{model}{A list containing information about the conformal prediction model.}
#' @author Xiaoqian Wang
#' @references Gibbs, I., and Candes, E. (2021). "Adaptive conformal inference under
#' distribution shift", \emph{Advances in Neural Information Processing Systems},
#' \bold{34}, 1660--1672.
#' @examples
#' 
#' library(forecast)
#' 
#' far2 <- function(x, h, level){
#'   Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level = level)
#' }
#' cvfc <- CVforecast(lynx, far2, h = 3, forward = TRUE, window = 30)
#' 
#' # ACP with symmetric conformity scores and rolling calibration sets
#' acpfc <- ACP(cvfc, symmetric = TRUE, gamma = 0.01,
#'              ncal = 30, rolling = TRUE)
#' 
#' @export
ACP <- function(object, alpha = 1 - 0.01 * object$level, gamma = 0.005,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                quantiletype = 1, na.rm = FALSE, ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (gamma < 0)
    stop("the step size parameter gamma should be positive")
  if (ncal < 10)
    stop("length of calibration period should at least be 10")
  if (!quantiletype %in% 1:9)
    stop("quantiletype is invalid. It must be in 1:9.")
  
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
  lower <- upper <- rep(list(namatrix), length(alpha))
  names(lower) <- names(upper) <- paste0(level, "%")
  alphat <- alphat_lower <- alphat_upper <- lower
  
  out <- list(
    x = object$x
  )
  
  for (h in seq(horizon)) {
    first_non_na <- (!is.na(errors[, h])) |> which() |> min()
    last_non_na <- (!is.na(errors[, h])) |> which() |> max()
    if (last_non_na < first_non_na + ncal - 1L)
      stop("errors in the input object is not long enough for calibration")
    indx <- seq(first_non_na + ncal - 1L, last_non_na + h - 1L, by = 1L)
    
    alphat_h <- alphat_lower_h <- alphat_upper_h <-
      errt_h <- errt_lower_h <- errt_upper_h <-
      matrix(NA_real_, nrow = NROW(errors), ncol = length(alpha))
    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, first_non_na, t - ncal + 1L),
        end = ifelse(t <= last_non_na, t, last_non_na))
      
      if (symmetric) {
        if (t == indx[1])
          alphat_h[t+1, ] <- alpha
        
        # Compute sample quantiles
        # (alpha_{t+1} is used to calculate sample quantiles from errors until t)
        q_lo <- q_up <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alphat_h[t+1, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        
        # Compute errt
        errt_h[t+1, ] <- abs(errors[t+1, h]) > q_lo
        outl <- which(alphat_h[t+1, ] >= 1)
        outs <- which(alphat_h[t+1, ] <= 0)
        errt_h[t+1, outl] <- TRUE
        errt_h[t+1, outs] <- FALSE
        
        if (t == tail(indx, 1))
          next
        if (t < last_non_na) {
          # Update alpha
          alphat_h[t+2, ] <- alphat_h[t+1, ] + gamma*(alpha - errt_h[t+1, ])
        } else {
          # Keep alpha unchanged after errt is not available
          alphat_h[t+2, ] <- alphat_h[t+1, ]
        }
      } else {
        if (t == indx[1])
          alphat_lower_h[t+1, ] <- alphat_upper_h[t+1, ] <- alpha/2
        q_lo <- ggdist::weighted_quantile(
          x = -c(errors_subset, Inf),
          probs = 1 - alphat_lower_h[t+1, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        q_up <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alphat_upper_h[t+1, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        
        # Compute errt
        errt_lower_h[t+1, ] <- (-errors[t+1, h]) > q_lo
        errt_lower_h[t+1, which(alphat_lower_h[t+1, ] >= 1)] <- TRUE
        errt_lower_h[t+1, which(alphat_lower_h[t+1, ] <= 0)] <- FALSE
        
        errt_upper_h[t+1, ] <- errors[t+1, h] > q_up
        errt_upper_h[t+1, which(alphat_upper_h[t+1, ] >= 1)] <- TRUE
        errt_upper_h[t+1, which(alphat_upper_h[t+1, ] <= 0)] <- FALSE
        
        if (t == tail(indx, 1))
          next
        if (t < last_non_na) {
          # Update alpha
          alphat_lower_h[t+2, ] <- alphat_lower_h[t+1, ] +
            gamma*(alpha/2 - errt_lower_h[t+1, ])
          alphat_upper_h[t+2, ] <- alphat_upper_h[t+1, ] +
            gamma*(alpha/2 - errt_upper_h[t+1, ])
        } else {
          # Keep alpha unchanged
          alphat_lower_h[t+2, ] <- alphat_lower_h[t+1, ]
          alphat_upper_h[t+2, ] <- alphat_upper_h[t+1, ]
        }
      }
      for (i in seq(length(alpha))) {
        lbl <- paste0(level[i], "%")
        lower[[lbl]][t+1, h] <- pf[t+1, h] - q_lo[i]
        upper[[lbl]][t+1, h] <- pf[t+1, h] + q_up[i]
      }
    }
    for (i in seq(length(alpha))) {
      alphat[[i]][, h] <- alphat_h[, i]
      alphat_lower[[i]][, h] <- alphat_lower_h[, i]
      alphat_upper[[i]][, h] <- alphat_upper_h[, i]
    }
  }
  if (h == 1) {
    lower <- lapply(lower, function(lo) lo[, 1L])
    upper <- lapply(lower, function(up) up[, 1L])
    alphat <- lapply(alphat, function(alp) alp[, 1L])
    alphat_lower <- lapply(alphat_lower, function(alp_lo) alp_lo[, 1L])
    alphat_upper <- lapply(alphat_upper, function(alp_up) alp_up[, 1L])
  }
  out$method <- paste("ACP")
  out$mean <- object$mean
  out$errors <- object$errors
  out$lower <- lower
  out$upper <- upper
  out$level <- level
  out$model$method <- out$method
  out$model$call <- match.call()
  out$model$alpha <- alpha
  out$model$gamma <- gamma
  out$model$symmetric <- symmetric
  out$model$ncal <- ncal
  out$model$rolling <- rolling
  out$model$quantiletype <- quantiletype
  if (symmetric) {
    out$model$alpha_update <- alphat
  } else {
    out$model$alpha_update <- list(lower = alphat_lower, upper = alphat_upper)
  }
  
  return(structure(out, class = "CPforecast"))
}

print.ACP <- function(y, ...) {
  x <- y$model
  cat(paste(x$method, "\n\n"))
  if (!is.null(x$call)) {
    cat(paste("Call:\n"))
    for (i in 1:length(deparse(x$call))) {
      cat(paste("", deparse(x$call)[i]), "\n")
    }
    cat(paste("\n"))
  }
  
  cat("  Conformal prediction settings:\n")
  cat(paste("    symmetric =", x$symmetric, "\n"))
  cat(paste("    rolling =", x$rolling, "\n"))
  cat(paste("   ", ifelse(rolling, "window =", "initial window ="), x$ncal, "\n"))
  
  cat("\n  Alpha update process:\n")
  cat(paste("    gamma =", x$gamma, "\n"))
  
  cat("\n  Sample quantiles:\n")
  cat(paste("    type =", x$quantiletype, "\n"))
  cat(paste("    weights = equal weights\n"))
}
