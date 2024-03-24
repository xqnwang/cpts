#' Adaptive conformal prediction method for time series forecasting
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
#' # Simulation series AR(2)
#' set.seed(0)
#' series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
#' series <- as.numeric(series)
#' 
#' # Setup
#' far2 <- function(x, h, level) {
#'   Arima(x, order = c(2, 0, 0)) |> 
#'     forecast(h = h, level)
#' }
#' 
#' # Cross-validation forecasting
#' fc <- CVforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
#'                  forward = TRUE, window = 100, initial = 1)
#' 
#' # ACP with symmetric conformity scores and rolling calibration sets
#' acpfc <- ACP(fc, symmetric = FALSE, gamma = 0.005, ncal = 100, rolling = TRUE)
#' 
#' @importFrom ggdist weighted_quantile
#' @export
acp <- function(object, alpha = 1 - 0.01 * object$level, gamma = 0.005,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                quantiletype = 1, na.rm = TRUE, ...) {
  # Check inputs
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
  pf <- ts(as.matrix(object$MEAN),
           start = start(object$MEAN),
           frequency = frequency(object$MEAN))
  errors <- ts(as.matrix(object$ERROR),
               start = start(object$ERROR),
               frequency = frequency(object$ERROR))
  horizon <- ncol(pf)
  n <- nrow(pf)
  
  if (ncal > nrow(errors))
    stop("`ncal` is larger than the number of rows in object$ERROR")
  
  namatrix <- ts(matrix(NA_real_, nrow = n, ncol = horizon),
                 start = start(pf),
                 frequency = frequency(pf))
  colnames(namatrix) <- paste0("h=", seq(horizon))
  lower <- upper <- alphat <- alphat_lower <- alphat_upper <-
    `names<-` (rep(list(namatrix), length(alpha)),
               paste0(level, "%"))
  
  out <- list(
    x = object$x,
    series = object$series
  )
  
  for (h in seq(horizon)) {
    indx <- seq(ncal, n - h, by = 1L)
    
    alphat_h <- alphat_lower_h <- alphat_upper_h <-
      errt_h <- errt_lower_h <- errt_upper_h <-
      matrix(NA_real_, nrow = n, ncol = length(alpha))
    
    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, 1, t - ncal + 1L),
        end = t)
      
      if (symmetric) {
        if (t == indx[1])
          alphat_h[t+h, ] <- alpha
        
        # Compute sample quantiles
        q_lo <- q_up <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alphat_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        
        # Compute errt
        errt_h[t+h, ] <- tryCatch(
          {abs(errors[t+h, h]) > q_lo},
          error = function(e) return(NA_real_)
        )
        outl <- which(alphat_h[t+h, ] >= 1)
        outs <- which(alphat_h[t+h, ] <= 0)
        errt_h[t+h, outl] <- TRUE
        errt_h[t+h, outs] <- FALSE
        
        if (t < tail(indx, 1)) {
          if (all(is.na(errt_h[t+1, ]))) {
            # Keep alpha unchanged
            alphat_h[t+h+1, ] <- alphat_h[t+h, ]
          } else {
            # Update alpha
            alphat_h[t+h+1, ] <- alphat_h[t+h, ] + gamma*(alpha - errt_h[t+1, ])
          }
        }
      } else {
        if (t == indx[1])
          alphat_lower_h[t+h, ] <- alphat_upper_h[t+h, ] <- alpha/2
        
        # Compute sample quantiles
        q_lo <- ggdist::weighted_quantile(
          x = -c(errors_subset, Inf),
          probs = 1 - alphat_lower_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        q_up <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alphat_upper_h[t+h, ],
          type = quantiletype,
          na.rm = na.rm,
          ...)
        
        # Compute errt
        errt_lower_h[t+h, ] <- tryCatch(
          {(-errors[t+h, h]) > q_lo},
          error = function(e) return(NA_real_)
        )
        errt_lower_h[t+h, which(alphat_lower_h[t+h, ] >= 1)] <- TRUE
        errt_lower_h[t+h, which(alphat_lower_h[t+h, ] <= 0)] <- FALSE
        
        errt_upper_h[t+h, ] <- tryCatch(
          {(errors[t+h, h]) > q_up},
          error = function(e) return(NA_real_)
        )
        errt_upper_h[t+h, which(alphat_upper_h[t+h, ] >= 1)] <- TRUE
        errt_upper_h[t+h, which(alphat_upper_h[t+h, ] <= 0)] <- FALSE
        
        if (t < tail(indx, 1)) {
          if (any(is.na(errt_lower_h[t+1, ])) || any(is.na(errt_upper_h[t+1, ]))) {
            # Keep alpha unchanged
            alphat_lower_h[t+h+1, ] <- alphat_lower_h[t+h, ]
            alphat_upper_h[t+h+1, ] <- alphat_upper_h[t+h, ]
          } else {
            # Update alpha
            alphat_lower_h[t+h+1, ] <- alphat_lower_h[t+h, ] + gamma*(alpha/2 - errt_lower_h[t+1, ])
            alphat_upper_h[t+h+1, ] <- alphat_upper_h[t+h, ] + gamma*(alpha/2 - errt_upper_h[t+1, ])
          }
        }
      }
      for (i in seq(length(alpha))) {
        lower[[i]][t+h, h] <- pf[t+h, h] - q_lo[i]
        upper[[i]][t+h, h] <- pf[t+h, h] + q_up[i]
      }
    }
    for (i in seq(length(alpha))) {
      alphat[[i]][, h] <- alphat_h[, i]
      alphat_lower[[i]][, h] <- alphat_lower_h[, i]
      alphat_upper[[i]][, h] <- alphat_upper_h[, i]
    }
  }
  
  out$method <- paste("acp")
  out$cp_times <- length(indx)
  out$MEAN <- object$MEAN
  out$ERROR <- object$ERROR
  out$LOWER <- lower
  out$UPPER <- upper
  out$level <- level
  out$call <- match.call()
  if ("mean" %in% names(object)) {
    out$mean <- object$mean
    out$lower <- extract_final(lower, nrow = n, ncol = horizon, bench = out$mean)
    out$upper <- extract_final(upper, nrow = n, ncol = horizon, bench = out$mean)
  }
  out$model$method <- out$method
  out$model$call <- match.call()
  out$model$alpha <- alpha
  out$model$gamma <- gamma
  out$model$symmetric <- symmetric
  if (symmetric) {
    out$model$alpha_update <- alphat
  } else {
    out$model$alpha_update <- list(lower = alphat_lower, upper = alphat_upper)
  }
  
  return(structure(out, class = c("acp", "cpforecast", "forecast")))
}
