#' @export
width <- function(object, ..., level = 95, includemedian = FALSE, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (level > 0 && level < 1) {
    level <- 100 * level
  } else if (level < 0 || level > 99.99) {
    stop("confidence limit out of range")
  }
  dots <- rlang::dots_list(...)
  if (missing(object)) {
    if (any(!(c("LOWER", "UPPER") %in% names(dots))))
      stop("LOWER, and UPPER are required for coverage calculation")
  } else {
    if (any(!(c("LOWER", "UPPER") %in% names(object))))
      stop("LOWER, and UPPER are required for coverage calculation")
    if (!(level %in% object$level))
      stop("no interval forecasts of target confidence level in object")
    levelname <- paste0(level, "%")
    LOWER <- object$LOWER[[levelname]]
    UPPER <- object$UPPER[[levelname]]
  }
  lower <- LOWER
  upper <- UPPER
  horizon <- ncol(lower)
  period <- frequency(lower)
  
  # Match time
  tspl <- tsp(lower)
  tspu <- tsp(upper)
  start <- max(tspl[1], tspu[1])
  end <- min(tspl[2], tspu[2])
  
  lower <- window(lower, start = start, end = end)
  upper <- window(upper, start = start, end = end)
  n <- nrow(lower)
  
  # Width matrix
  widmat <- (upper- lower) |>
    ts(start = start, end = end, frequency = period)
  colnames(widmat) <- colnames(lower)
  
  out <- list(
    width = widmat
  )
  
  # Mean coverage
  out$mean <- apply(widmat, 2, mean, na.rm = na.rm)
  
  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    out$rollmean <- apply(widmat, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
  }
  
  # Median
  if (includemedian) {
    # Mean median
    out$median <- apply(widmat, 2, median, na.rm = na.rm)
    
    # Rolling median coverage
    if (!is.null(window)) {
      if (window >= n)
        stop("the `window` argument should be smaller than the total period of interest")
      out$rollmedian <- apply(widmat, 2, zoo::rollmedian, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
    }
  }
  
  return(structure(out, class = "width"))
}

print.width <- function(x, ...) {
  cat("Mean width:\n")
  print(x$mean)
  
  if ("median" %in% names(x)) {
    cat("\nMedian width:\n")
    print(x$median)
  }
}
