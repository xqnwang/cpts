#' @importFrom stats window
#' @import zoo rollmean
#' @export
coverage <- function(object, ..., level = 95, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (level > 0 && level < 1) {
    level <- 100 * level
  } else if (level < 0 || level > 99.99) {
    stop("confidence limit out of range")
  }
  dots <- rlang::dots_list(...)
  if (missing(object)) {
    if (any(!(c("x", "LOWER", "UPPER") %in% names(dots))))
      stop("x, LOWER, and UPPER are required for coverage calculation")
  } else {
    if (any(!(c("x", "LOWER", "UPPER") %in% names(object))))
      stop("x, LOWER, and UPPER are required for coverage calculation")
    if (!(level %in% object$level))
      stop("no interval forecasts of target confidence level in object")
    levelname <- paste0(level, "%")
    x <- object$x
    LOWER <- object$LOWER[[levelname]]
    UPPER <- object$UPPER[[levelname]]
  }
  lower <- LOWER
  upper <- UPPER
  horizon <- ncol(lower)
  period <- frequency(object$x)
  x <- ts(matrix(rep(object$x, horizon), ncol = horizon, byrow = FALSE),
          start = start(object$x),
          frequency = period)
  
  # Match time
  tspx <- tsp(x)
  tspl <- tsp(lower)
  tspu <- tsp(upper)
  start <- max(tspx[1], tspl[1], tspu[1])
  end <- min(tspx[2], tspl[2], tspu[2])
  
  x <- window(x, start = start, end = end)
  lower <- window(lower, start = start, end = end)
  upper <- window(upper, start = start, end = end)
  n <- nrow(x)
  
  # If coverage matrix
  covmat <- (lower <= x & x <= upper) |>
    ts(start = start, end = end, frequency = period)
  colnames(covmat) <- colnames(lower)
  
  # Mean coverage
  covmean <- apply(covmat, 2, mean, na.rm = na.rm)
  
  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    covrmean <- apply(covmat, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
      ts(end = end, frequency = period)
  }
  
  out <- list(
    mean = covmean,
    ifinn = covmat
  )
  if (!is.null(window)) out <- append(out, list(rolling = covrmean))
  return(structure(out, class = "coverage"))
}

print.coverage <- function(x, ...) {
  print(x$mean)
}
