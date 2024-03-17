#' Time series cross-validation forecasting
#' 
#' \code{CVforecast} computes forecasts and other information obtained by applying
#' \code{forecastfun} to subsets of the time series \code{y} using a
#' rolling forecast origin.
#' 
#' Let \code{y} contain the time series \eqn{y_1,\dots,y_T}{y[1:T]}. Let initial
#' period be \eqn{t_0} and \code{forward = TRUE}. If window is NULL,
#' \code{forecastfun} is applied successively to the time series
#' \eqn{y_{1},\dots,y_t}{y[1:t]}, for \eqn{t=t_0,\dots,T}, making forecasts
#' \eqn{\hat{y}_{t+h|t}}{f[t+h]}. If window is not NULL and of length \eqn{t_w},
#' \code{forecastfun} is applied successively to the time series
#' \eqn{y_{t-t_w+1},\dots,y_{t}}{y[(t-t_w+1):t)]}, for \eqn{t=\max(t_0, t_w),\dots,T}.
#' If \code{forward = FALSE}, the last observation used for training is \eqn{y_{T-1}}.
#' 
#' @param y Univariate time series.
#' @param forecastfun Function to return an object of class \code{forecast}.
#' Its first argument must be a univariate time series, and it must have an
#' argument \code{h} for the forecast horizon and an argument \code{level} for
#' the confidence level for prediction intervals. If exogenous predictors are used,
#' then it must also have \code{xreg} and \code{newxreg} arguments corresponding
#' to the training and test periods, respectively.
#' @param h Forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param forward If \code{TRUE}, the last observation used for forecasting is \eqn{y_T}.
#' Otherwise, the last observation used for forecasting is \eqn{y_{T-1}}.
#' @param window Length of the rolling window. If NULL, a rolling window will not be used.
#' @param xreg Exogenous predictor variables passed to \code{forecastfun} if required.
#' @param newxreg Exogenous predictor variables passed to \code{forecastfun} if
#' required when \code{forward = TRUE}.
#' @param initial Initial period of the time series where no cross-validation is performed.
#' @param ... Other arguments are passed to \code{forecastfun}.
#' @return An object of class "\code{CVforecast}".
#' 
#' An object of class "\code{CVforecast}" is a list containing the following elements:
#' \item{x}{The original time series.}
#' \item{method}{The name of the forecasting method as a character string.}
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
#' \item{model}{A list containing information about the latest fitted model.}
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
CVforecast <- function(y, forecastfun, h = 1, level = c(80, 95),
                       forward = TRUE, window = NULL, xreg = NULL, newxreg = NULL,
                       initial = 1, ...) {
  y <- as.ts(y)
  n <- length(y)
  level <- sort(level)
  
  if (h <= 0) 
    stop("forecast horizon out of bounds")
  if (initial <= 0)
    stop("initial period should be at least of length 1")
  if (initial >= n)
    stop("initial period too long")
  if (!is.null(window)) {
    if (window >= n)
      stop("window period too long")
  }
  if (forward) {
    N <- n + h
    nlast <- n
  } else {
    N <- n + h - 1L
    nlast <- n - 1L
  }
  
  if (!is.null(xreg)) {
    # Make xreg a ts object to allow easy subsetting later
    xreg <- ts(as.matrix(xreg))
    if (NROW(xreg) != length(y))
      stop("xreg must be of the same size as y")
    
    if (forward && !is.null(newxreg)) {
      newxreg <- ts(as.matrix(newxreg))
      if (NROW(newxreg) < h)
        stop("newxreg must be of the same size as h")
      if (NROW(newxreg) > h) {
        message("only first h rows of newxreg are being used for forecasting")
        newxreg <- ts(as.matrix(newxreg[seq(h),]))
      }
      # Pad xreg with newxreg
      xreg <- ts(rbind(xreg, newxreg),
                 start = start(y),
                 frequency = frequency(y))
    } else {
      # Pad xreg with NAs
      xreg <- ts(rbind(xreg, matrix(NA_real_, nrow = h, ncol = NCOL(xreg))),
                 start = start(y),
                 frequency = frequency(y))
    }
  }
  
  if (is.null(window)) {
    indx <- seq(initial, nlast, by = 1L)
  } else {
    indx <- seq(max(window, initial), nlast, by = 1L)
  }
  
  pf <- err <- ts(matrix(NA_real_, nrow = N, ncol = h), 
                  start = start(y), 
                  frequency = frequency(y))
  colnames(pf) <- colnames(err) <- paste0("h=", seq(h))
  lower <- rep(list(pf), length(level))
  upper <- rep(list(pf), length(level))
  names(lower) <- names(upper) <- paste0(level, "%")
  out <- list(
    x = y
  )
  
  for (i in indx) {
    y_subset <- subset(
      y,
      start = ifelse(is.null(window), 1L, i - window + 1L),
      end = i)
    
    if (is.null(xreg)) {
      fc <- try(suppressWarnings(
        forecastfun(y_subset, h = h, level = level, ...)
        ), silent = TRUE)
    } else {
      xreg_subset <- subset(
        xreg,
        start = ifelse(is.null(window), 1L, i - window + 1L),
        end = i)
      xreg_future <- subset(
        xreg,
        start = i + 1L,
        end = i + h)
      fc <- try(suppressWarnings(
        forecastfun(y_subset, h = h, level = level,
                    xreg = xreg_subset, newxreg = xreg_future, ...)
        ), silent = TRUE)
    }
    
    if (!is.element("try-error", class(fc))) {
      pf[i,] <- fc$mean
      err[i,] <- y[i + seq(h)] - fc$mean
      for (l in level) {
        lower[[paste0(l, "%")]][i,] <- fc$lower[seq(h), paste0(l, "%")]
        upper[[paste0(l, "%")]][i,] <- fc$upper[seq(h), paste0(l, "%")]
      }
    }
  }
  
  out$method <- as.character(fc$model$call)[1]
  out$mean <- lagmatrix(pf, seq(h)) |> subset(start = indx[1] + 1L)
  out$errors <- lagmatrix(err, seq(h)) |> subset(start = indx[1] + 1L)
  out$lower <- lapply(lower,
                      function(low) lagmatrix(low, seq(h)) |>
                        subset(start = indx[1] + 1L))
  out$upper <- lapply(upper,
                      function(up) lagmatrix(up, seq(h)) |>
                        subset(start = indx[1] + 1L))
  out$level <- level
  # The final forecasting model output in the for loop
  out$model <- fc$model
  
  return(structure(out, class = "CVforecast"))
}

lagmatrix <- function(x, k) {
  # Ensure 'x' is a matrix
  if (!is.matrix(x))
    stop("ensure x is a matrix")
  if (NCOL(x) != length(k))
    stop("k must have the same number of columns as x")
  
  n <- NROW(x)
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