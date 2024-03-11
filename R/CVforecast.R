#' Time series cross-validation forecasting
#' 
#' \code{CVforecast} computes forecasts and other information obtained by applying
#' \code{forecastfunction} to subsets of the time series \code{y} using a
#' rolling forecast origin.
#' 
#' Let \code{y} contain the time series \eqn{y_1,\dots,y_T}{y[1:T]}. Let initial
#' period be \eqn{t_0} and \code{forward = TRUE}. If window is NULL,
#' \code{forecastfunction} is applied successively to the time series
#' \eqn{y_{1},\dots,y_t}{y[1:t]}, for \eqn{t=t_0+1,\dots,T}, making predictions
#' \eqn{\hat{y}_{t+h|t}}{f[t+h]}. If window is not NULL and of length \eqn{t_w},
#' \code{forecastfunction} is applied successively to the time series
#' \eqn{y_{t-t_w+1},\dots,y_{t}}{y[(t-t_w+1):t)]}, for \eqn{t=\max(t_0, t_w)+1,\dots,T}.
#' If \code{forward = FALSE}, the last observation used for training is \eqn{y_{T-1}}.
#' 
#' @param y Univariate time series
#' @param forecastfunction Function to return an object of class \code{forecast}.
#' Its first argument must be a univariate time series, and it must have an
#' argument \code{h} for the forecast horizon and an argument \code{level} for
#' the confidence level for prediction intervals. If exogenous predictors are used,
#' then it must also have \code{xreg} and \code{newxreg} arguments corresponding
#' to the training and test periods.
#' @param h Forecast horizon
#' @param level Confidence level for prediction intervals
#' @param PI 	If \code{TRUE}, prediction intervals are produced, otherwise only point
#' forecasts are calculated.
#' @param forward If \code{TRUE}, the last observation used for forecasting is \eqn{y_T}.
#' Otherwise, the last observation used for forecasting is \eqn{y_{T-1}}.
#' @param window Length of the rolling window, if NULL, a rolling window will not be used.
#' @param xreg Exogeneous predictor variables passed to the forecast function if required
#' @param newxreg Exogeneous predictor variables passed to the forecast function if required when \code{forward = TRUE}
#' @param initial Initial period of the time series where no cross-validation is performed
#' @param ... Other arguments are passed to \code{forecastfunction}
#' @return An object of class "\code{CVforecast}".
#' 
#' An object of class \code{"CVforecast"} is a list containing at least the
#' following elements:
#' \item{model}{A list containing information about the latest fitted model}
#' \item{method}{The name of the latest forecasting method as a character string}
#' \item{mean}{Point forecasts as a time series for \eqn{h=1}. For \eqn{h>1},
#' they are returned as a time series matrix with the \eqn{h}th column
#' containing point forecasts for forecast horizon \eqn{h}. The time index
#' corresponds to the forecast period.}
#' \item{errors}{Forecast errors given by
#' \eqn{e_{t+h} = y_{t+h}-\hat{y}_{t+h|t}}{e[t+h] = y[t+h]-f[t+h]}}
#' \item{lower}{Lower limits for prediction intervals}
#' \item{upper}{Upper limits for prediction intervals}
#' \item{level}{The confidence values associated with the prediction intervals}
#' \item{x}{The original time series}
#' \item{residuals}{Residuals from the latest fitted model.
#' For models with additive errors, the residuals are x - fitted values. For
#' models with multiplicative errors, the residuals are equal to x /(fitted
#' values) - 1.}
#' \item{fitted}{Fitted values from the latest fitted model (one-step forecasts)}
#' @author Xiaoqian Wang
#' @examples
#' 
#' library(forecast)
#' 
#' # Fit an AR(2) model to each rolling origin subset until forecast origin T-1
#' far2 <- function(x, h, level){
#'   Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level = level)
#' }
#' cvfc <- CVforecast(lynx, far2, h = 1, forward = FALSE)
#' 
#' # Fit the same model with a rolling forecast origin until T
#' cvfc <- CVforecast(lynx, far2, h = 1, forward = TRUE)
#' 
#' # Fit the same model with a rolling window of length 30, a rolling forecast origin until T, and a forecast horizon set to 3
#' cvfc <- CVforecast(lynx, far2, h = 3, forward = TRUE, window = 30)
#' 
#' # Example with exogenous predictors and a rolling window of length 30
#' far2_xreg <- function(x, h, level, xreg, newxreg) {
#'   Arima(x, order = c(2, 0, 0), xreg = xreg) |> 
#'     forecast(xreg = newxreg, level)
#' }
#' 
#' y <- ts(rnorm(100))
#' xreg <- matrix(rnorm(200), ncol = 2)
#' newxreg <- matrix(rnorm(10), ncol = 2)
#' cvfc <- CVforecast(y, far2_xreg, h = 5, forward = TRUE, window = 30,
#'                    xreg = xreg, newxreg = newxreg)
#' 
#' cvfc <- CVforecast(y, far2_xreg, h = 5, forward = FALSE, window = 30,
#'                    xreg = xreg, newxreg = newxreg)
#' 
#' @export
CVforecast <- function(y, forecastfunction, h = 1, level = c(80, 95), PI = TRUE,
                       forward = TRUE, window = NULL, xreg = NULL, newxreg = NULL,
                       initial = 0, ...) {
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
    ntr <- n
  } else {
    N <- n + h - 1L
    ntr <- n - 1L
  }
  
  if (!is.null(xreg)) {
    # Make xreg a ts object to allow easy subsetting later
    xreg <- ts(as.matrix(xreg))
    if (NROW(xreg) != length(y))
      stop("xreg must be of the same size as y")
    
    if (forward && !is.null(newxreg)) {
      newxreg <- ts(as.matrix(newxreg))
      if (NROW(newxreg) != h)
        stop("newxreg must be of the same size as h")
      # Pad xreg with newxreg
      xreg <- ts(rbind(xreg, newxreg),
                 start = start(y),
                 frequency = frequency(y))
    } else {
      # Pad xreg with NAs
      xreg <- ts(rbind(xreg, matrix(NA, nrow = h, ncol = NCOL(xreg))),
                 start = start(y),
                 frequency = frequency(y))
    }
  }
  
  if (is.null(window)) {
    indx <- seq(1 + initial, ntr, by = 1L)
  } else {
    indx <- seq(1 + max(window, initial), ntr, by = 1L)
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
        forecastfunction(y_subset, h = h, level = level, ...)
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
        forecastfunction(y_subset, h = h, level = level,
                         xreg = xreg_subset, newxreg = xreg_future, ...)
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
  # The final forecasting model output in the for loop
  out$method <- fc$method
  out$model <- fc$model
  out$fitted <- fc$fitted
  out$series <- fc$series
  out$residuals <- fc$residuals
  
  if (h == 1) {
    out$mean <- out$mean[, 1L]
    out$errors <- out$errors[, 1L]
    out$lower <- lapply(out$lower, function(lo) lo[, 1L])
    out$upper <- lapply(upper$lower, function(up) up[, 1L])
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

