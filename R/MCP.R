#' Multistep-ahead conformal prediction
#' 
#' \code{MCP} computes prediction intervals and other information obtained by
#' applying the multistep-ahead conformal prediction method. The method can only
#' deal with asymmetric conformity scores, i.e., forecast errors.
#' 
#' @param object An object of class "\code{CVforecast}". It must have an argument
#' \code{x} for original univariate time series, an argument \code{mean} for
#' point forecasts and \code{errors} for forecast errors. See the results of a call
#' to \code{\link{CVforecast}}.
#' @param alpha A numeric vector of target levels \eqn{1-\alpha}.
#' @param ncal Length of the calibration set. If \code{rolling=FALSE}, it denotes
#' the initial period of calibration sets. If \code{rolling=TRUE}, it indicates
#' the period of each rolling calibration set.
#' @param rolling If \code{TRUE}, a rolling window strategy will be used to
#' generate the calibration set for learning rate update and scorecasting.
#' Otherwise, expanding window will be used.
#' @param integrate If \code{TRUE}, error integration will be included in the
#' update process.
#' @param scorecast If \code{TRUE}, scorecasting will be included in the update
#' process as a forecast combination of a MA(h-1) model and a linear regression model.
#' @param lr Initial learning rate used for quantile tracking.
#' @param Tg The time set to achieve the target absolute guarantee before this.
#' @param delta The target absolute guarantee is \eqn{1-\alpha-\delta} coverage.
#' @param Csat A positive constant ensuring that by time \code{Tg}, an absolute
#' guarantee is of at least \eqn{1-\alpha-\delta} coverage.
#' @param KI A positive constant to place the integrator on the same scale as the scores.
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
#' # PID setup
#' Tg <- 1000
#' delta <- 0.01
#' lr <- 0.1
#' KI <- 2
#' Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
#' 
#' # MCP with integrator and scorecaster
#' mcpfc <- MCP(fc, ncal = 100, rolling = TRUE,
#'              integrate = TRUE, scorecast = TRUE,
#'              lr = lr, KI = KI, Csat = Csat)
#' 
#' @importFrom stats lm
#' @importFrom forecast meanf Arima forecast
#' @export
mcp <- function(object, alpha = 1 - 0.01 * object$level,
                ncal = 10, rolling = FALSE,
                integrate = TRUE, scorecast = TRUE,
                lr = 0.1, Tg = NULL, delta = NULL,
                Csat = 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg)),
                KI = abs(object$errors) |> max(na.rm = TRUE), ...) {
  # Check inputs
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (ncal < 10)
    stop("Length of calibration period should at least be 10")
  if (!is.null(Tg) && !is.null(delta)) {
    if (!is.null(Csat))
      warning("Csat is replaced by calculation using Tg and delta")
    Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
  }
  
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
  
  namatrix <- `colnames<-` (
    ts(matrix(NA_real_, nrow = n, ncol = horizon),
       start = start(pf),
       frequency = frequency(pf)),
    paste0("h=", seq(horizon)))
  nalist <- `names<-` (
    rep(list(namatrix), length(alpha)),
    paste0(level, "%"))
  
  lower <- upper <- nalist
  lrmat <- namatrix
  if (integrate)
    integrator_lower <- integrator_upper <- nalist
  if (scorecast)
    scorecaster_lower <- scorecaster_upper <- namatrix
  
  out <- list(
    x = object$x,
    series = object$series
  )
  
  for (h in seq(horizon)) {
    indx <- seq(h, nrow(errors), by = 1L)
    
    errt_lower_h <- errt_upper_h <-
      integ_lower_h <- integ_upper_h <-
      q_lo_h <- q_up_h <-
      matrix(NA_real_, nrow = n, ncol = length(alpha))
    qts_lower_h <- qts_upper_h <-
      qs_lower_h <- qs_upper_h <- matrix(0, nrow = n, ncol = length(alpha))
    for (t in indx) {
      t_burnin <- max(t - ncal + 1L, h)
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, h, t_burnin),
        end = t)
      
      # Calculate errt
      errt_lower_h[t, ] <- (-errors[t, h]) > qs_lower_h[t, ]
      errt_upper_h[t, ] <- errors[t, h] > qs_upper_h[t, ]
      
      # Learning rate (same for the upper and lower bounds)
      lrmat[t, h] <- ifelse(
        length(errors_subset) <= 1,
        lr,
        lr*(max(errors_subset) - min(errors_subset)))
      
      # Update quantile tracking
      qts_lower_h[t+h, ] <- qts_lower_h[t+h-1, ] + lrmat[t, h] * (errt_lower_h[t, ] - alpha/2)
      qts_upper_h[t+h, ] <- qts_upper_h[t+h-1, ] + lrmat[t, h] * (errt_upper_h[t, ] - alpha/2)
      
      # Update integrator
      if (integrate) {
        el <- errt_lower_h[h:t, ] |>
          matrix(ncol = length(alpha))
        integrator_lower_arg <- apply(el, 2, sum) - nrow(el)*alpha/2
        integ_lower_h[t+h, ] <- sapply(
          1:length(alpha),
          function(i) ifelse(
            nrow(el) == 1,
            0,
            saturation_fn_log(integrator_lower_arg[i], nrow(el), Csat, KI)))
        
        eu <- errt_upper_h[h:t, ] |>
          matrix(ncol = length(alpha))
        integrator_upper_arg <- apply(eu, 2, sum) - nrow(eu)*alpha/2
        integ_upper_h[t+h, ] <- sapply(
          1:length(alpha),
          function(i) ifelse(
            nrow(eu) == 1,
            0,
            saturation_fn_log(integrator_upper_arg[i], nrow(eu), Csat, KI)))
      }
      
      # Update scorecaster
      do_scorecast <- (scorecast && t >= (ncal+h-1))
      if (do_scorecast) {
        if (h == 1) {
          model <- forecast::meanf(errors_subset, h = h)
          scorecaster_lower[t+h, h] <- -as.numeric(model$mean)
          scorecaster_upper[t+h, h] <- as.numeric(model$mean)
        } else {
          model_MA <- forecast::Arima(errors_subset, order = c(0, 0, h-1)) |>
            forecast::forecast(h = h)
          model_LR <- lm(
            as.formula(paste0("V", h, " ~ .")), 
            data = setNames(
              as.data.frame(sapply(1:h, function(j) {
                subset(
                  errors[, j],
                  start = ifelse(!rolling, j, t-ncal+1-h+j),
                  end = t-h+j)
              })),
              paste0("V", 1:h))
          ) |>
            forecast::forecast(
              newdata = setNames(
                data.frame(sapply(1:(h-1),
                                  function(j) scorecaster_upper[t-h+1+j, j]) |> matrix(nrow = 1)),
                paste0("V", 1:(h-1))))
          scorecaster_lower[t+h, h] <- -(as.numeric(model_MA$mean[h]) + as.numeric(model_LR$mean))/2
          scorecaster_upper[t+h, h] <- (as.numeric(model_MA$mean[h]) + as.numeric(model_LR$mean))/2
        }
      }
      
      # Update the next quantile
      qs_lower_h[t+h, ] <- qts_lower_h[t+h, ] +
        ifelse(rep(integrate, length(alpha)), integ_lower_h[t+h, ], rep(0, length(alpha))) +
        rep(ifelse(do_scorecast,
                   ifelse(is.na(scorecaster_lower[t+h, h]), 0, scorecaster_lower[t+h, h]),
                   0),
            length(alpha))
      qs_upper_h[t+h, ] <- qts_upper_h[t+h, ] +
        ifelse(rep(integrate, length(alpha)), integ_upper_h[t+h, ], rep(0, length(alpha))) +
        rep(ifelse(do_scorecast,
                   ifelse(is.na(scorecaster_upper[t+h, h]), 0, scorecaster_upper[t+h, h]),
                   0),
            length(alpha))
      
      # PIs
      if (t >= (ncal+h-1)) {
        for (i in seq(length(alpha))) {
          lower[[i]][t+h, h] <- pf[t+h, h] - qs_lower_h[t+h, i]
          upper[[i]][t+h, h] <- pf[t+h, h] + qs_upper_h[t+h, i]
        }
      }
    }
    if (integrate) {
      for (i in seq(length(alpha))) {
        integrator_lower[[i]][, h] <- integ_lower_h[, i]
        integrator_upper[[i]][, h] <- integ_upper_h[, i]
      }
    }
  }
  
  out$method <- paste("mcp")
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
  out$model$symmetric <- symmetric
  out$model$integrate <- integrate
  out$model$scorecast <- scorecast
  out$model$lr <- lr
  out$model$Csat <- Csat
  out$model$KI <- KI
  out$model$lr_update <- lrmat
  if (integrate)
    out$model$integrator <- list(lower = integrator_lower, upper = integrator_upper)
  if (scorecast)
    out$model$scorecaster <- list(lower = scorecaster_lower, upper = scorecaster_upper)
  
  return(structure(out, class = c("mcp", "cpforecast", "forecast")))
}
