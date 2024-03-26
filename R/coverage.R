#' @importFrom stats window
#' @import zoo rollmean
#' @export
coverage <- function(object, level = object$level, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (any(!(c("x", "LOWER", "UPPER") %in% names(object))))
    stop("x, LOWER, and UPPER are required for coverage calculation")
  if (!(is.list(object$LOWER) && is.list(object$UPPER)))
    stop("LOWER and UPPER should be a list")
  if (all(level > 0 & level < 1)) {
    level <- 100 * level
  } else if (any(level < 0 | level > 99.99)) {
    stop("confidence limit out of range")
  }
  
  # Extract information of interest
  level <- sort(level)
  levelname <- paste0(level, "%")
  lower <- object$LOWER[levelname]
  upper <- object$UPPER[levelname]
  horizon <- ncol(lower[[1]])
  
  period <- frequency(object$x)
  x <- ts(matrix(rep(object$x, horizon), ncol = horizon, byrow = FALSE),
          start = start(object$x),
          frequency = period)
  
  # Match time
  tspx <- tsp(x)
  tspl <- tsp(lower[[1]])
  tspu <- tsp(upper[[1]])
  start <- max(tspx[1], tspl[1], tspu[1])
  end <- min(tspx[2], tspl[2], tspu[2])
  
  x <- window(x, start = start, end = end)
  lower <- lapply(lower, function(lo) window(lo, start = start, end = end))
  upper <- lapply(upper, function(up) window(up, start = start, end = end))
  n <- nrow(x)
  
  # If coverage matrix
  covmat <- `names<-` (lapply(levelname, function(i) {
    lo <- lower[[i]]
    up <- upper[[i]]
    xx <- x
    cov <- (lo <= xx & xx <= up)
    cov <- ts(cov, start = start, end = end, frequency = period)
    colnames(cov) <- colnames(lo)
    return(cov)
  }), levelname)
  
  # Mean coverage
  covmean <- sapply(covmat, function(cov) {
    apply(cov, 2, mean, na.rm = na.rm)
  })
  
  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    covrmean <- lapply(covmat, function(cov) {
      apply(cov, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
    })
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
