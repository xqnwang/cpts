#' @export
width <- function(object, ..., level = 95, window = NULL, na.rm = TRUE) {
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
  
  # Mean coverage
  widmean <- apply(widmat, 2, mean, na.rm = na.rm)
  
  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    widrmean <- apply(widmat, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
  }
  
  out <- list(
    mean = widmean,
    width = widmat
  )
  if (!is.null(window)) out <- append(out, list(rolling = widrmean))
  return(structure(out, class = "width"))
}

print.width <- function(x, ...) {
  print(x$mean)
}
