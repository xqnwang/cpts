#' Classical split conformal prediction methods for time series forecasting
#' 
#' \code{SCP} computes prediction intervals and other information obtained by
#' applying the classical split conformal prediction method.
#' 
#' Let \code{s} contain the \eqn{h}-step-ahead conformity scores
#' \eqn{s_{1+h|1},s_{2+h|2},\dots,s_{T|T-h}}. Let initial period \code{ncal} be
#' \eqn{t_c}.
#' 
#' If \code{symmetric=TRUE}, \eqn{s_{t+h|t}=|e_{t+h|t}|}. If \code{rolling=FALSE},
#' the \code{\link[ggdist]{weighted_quantile}} function is applied successively
#' to the calibration set \eqn{s_{1+h|1},\dots,s_{t|t-h}} to generate its
#' \eqn{\alpha}-quantile \eqn{\hat{q}_{t+h|t}}, for \eqn{t=t_c+h,\dots,T},
#' getting the prediction intervals \eqn{\hat{y}_{t+h|t} \pm \hat{q}_{t+h|t}}.
#' If \code{rolling=TRUE}, the \code{\link[ggdist]{weighted_quantile}} function
#' is applied successively to the calibration set
#' \eqn{s_{t-t_c+1|t-h-t_c+1},\dots,s_{t|t-h}} to generate \eqn{\hat{q}_{t+h|t}},
#' where \eqn{t=t_c+h,\dots,T}.
#' 
#' If \code{symmetric=FALSE}, \eqn{s_{t+h|t}^{u}=e_{t+h|t}} for upper limit, and
#' \eqn{s_{t+h|t}^{l} = -e_{t+h|t}} for lower limit. The \eqn{\alpha/2}-quantiles,
#' \eqn{\hat{q}_{t+h|t}^{u}} and \eqn{\hat{q}_{t+h|t}^{l}}, are generated for
#' upper limit and lower limit based on expanding/rolling calibration sets,
#' respectively. The final prediction interval is
#' \eqn{[\hat{y}_{t+h|t}-\hat{q}_{t+h|t}^{l}, \hat{y}_{t+h|t}+\hat{q}_{t+h|t}^{u}]}.
#' 
#' @param object An object of class "\code{CVforecast}". It must have an argument
#' \code{x} for original univariate time series, an argument \code{mean} for
#' point forecasts and \code{errors} for forecast errors. See the results of a call
#' to \code{\link{CVforecast}}.
#' @param alpha A numeric vector of target levels \eqn{1-\alpha}.
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
#' @param weightfun Function to return a vector of weights used for sample quantile
#' computation. Its first argument must be an integer indicating the number of
#' observations for which weights are generated. If \code{NULL}, equal weights
#' are used for sample quantile computation.
#' @param kess If \code{TRUE}, Kish's effective sample size is used for sample
#' quantile computation.
#' @param na.rm If \code{TRUE}, \code{NA} values are removed when calculating
#' sample quantile.
#' @param ... Other arguments are passed to \code{weightfun}.
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
#' far2 <- function(x, h, level){
#'   Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level = level)
#' }
#' cvfc <- CVforecast(lynx, far2, h = 3, forward = TRUE, window = 30)
#' 
#' # SCP with symmetric conformity scores and rolling calibration sets
#' scpfc <- SCP(cvfc, symmetric = TRUE,
#'              ncal = 30, rolling = TRUE)
#' 
#' # SCP with exponential weights for sample quantiles
#' expweight <- function(n, base = 0.99) c(base^((n-1):1), 1)
#' scpfc_exp <- SCP(cvfc, symmetric = TRUE,
#'                  ncal = 30, rolling = TRUE, weightfun = expweight)
#' 
#' @export
SCP <- function(object, alpha = 1 - 0.01 * object$level,
                symmetric = FALSE, ncal = 10, rolling = FALSE,
                quantiletype = 1, weightfun = NULL, kess = FALSE, na.rm = FALSE,
                ...) {
  if (any(alpha >= 1 | alpha <= 0))
    stop("alpha should be in (0, 1)")
  if (ncal < 10)
    stop("length of calibration period should at least be 10")
  if (!quantiletype %in% 1:9)
    stop("quantiletype is invalid. It must be in 1:9.")
  if (is.null(weightfun)) {
    weightfun <- function(n) rep(1, n + 1L)
  }
  # Kish's effective sample size for sample quantile computation
  if (kess) {
    kess <- function(w) sum(w)^2 / sum(w^2)
  } else {
    kess <- NULL
  }
  
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
  
  out <- list(
    x = object$x
  )
  
  for (h in seq(horizon)) {
    first_non_na <- (!is.na(errors[, h])) |> which() |> min()
    last_non_na <- (!is.na(errors[, h])) |> which() |> max()
    if (last_non_na < first_non_na + ncal - 1L)
      stop("errors in the input object is not long enough for calibration")
    indx <- seq(first_non_na + ncal - 1L, last_non_na + h - 1L, by = 1L)
    
    for (t in indx) {
      errors_subset <- subset(
        errors[, h],
        start = ifelse(!rolling, first_non_na, t - ncal + 1L),
        end = ifelse(t <= last_non_na, t, last_non_na))
      
      weight_subset <- weightfun(length(errors_subset) + 1L, ...)
      
      if (symmetric) {
        q_lo <- q_up <- ggdist::weighted_quantile(
          x = abs(c(errors_subset, Inf)),
          probs = 1 - alpha,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
      } else {
        q_lo <- ggdist::weighted_quantile(
          x = - c(errors_subset, Inf),
          probs = 1 - alpha/2,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
        q_up <- ggdist::weighted_quantile(
          x = c(errors_subset, Inf),
          probs = 1 - alpha/2,
          weights = weight_subset,
          n = kess,
          type = quantiletype,
          na.rm = na.rm)
      }
      for (i in seq(length(alpha))) {
        lbl <- paste0(level[i], "%")
        lower[[lbl]][t+1, h] <- pf[t+1, h] - q_lo[i]
        upper[[lbl]][t+1, h] <- pf[t+1, h] + q_up[i]
      }
    }
  }
  if (h == 1) {
    lower <- lapply(lower, function(lo) lo[, 1L])
    upper <- lapply(upper, function(up) up[, 1L])
  }
  out$method <- paste("SCP")
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
  out$model$quantiletype <- quantiletype
  out$model$weight <- ifelse(is.null(out$model$call$weightfun),
                             "Equal weights",
                             "User specified weights")
  out$model$kess <- out$model$call$kess
  
  return(structure(out, class = "CPforecast"))
}

print.SCP <- function(y, ...) {
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
  
  cat("\n  Sample quantiles:\n")
  cat(paste("    type =", x$quantiletype, "\n"))
  cat(paste("    weights =", tolower(x$weight), "\n"))
  cat(paste("    kess =", x$kess, "\n"))
}
