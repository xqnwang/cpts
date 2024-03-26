#' @export
width <- function(object, level = object$level, window = NULL, na.rm = TRUE) {
  # Check inputs
  if (any(!(c("LOWER", "UPPER") %in% names(object))))
    stop("LOWER, and UPPER are required for coverage calculation")
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
  period <- frequency(lower[[1]])
  
  # Match time
  tspl <- tsp(lower[[1]])
  tspu <- tsp(upper[[1]])
  start <- max(tspl[1], tspu[1])
  end <- min(tspl[2], tspu[2])
  
  lower <- lapply(lower, function(lo) window(lo, start = start, end = end))
  upper <- lapply(upper, function(up) window(up, start = start, end = end))
  n <- nrow(lower[[1]])
  
  # Width matrix
  widmat <- `names<-` (lapply(levelname, function(i) {
    wid <- upper[[i]] - lower[[i]]
    wid <- ts(wid, start = start, end = end, frequency = period)
    colnames(wid) <- colnames(upper[[i]])
    return(wid)
  }), levelname)
  
  # Mean coverage
  widmean <- sapply(widmat, function(wid) {
    apply(wid, 2, mean, na.rm = na.rm)
  })
  
  # Rolling mean coverage
  if (!is.null(window)) {
    if (window >= n)
      stop("the `window` argument should be smaller than the total period of interest")
    widrmean <- lapply(widmat, function(wid) {
      apply(wid, 2, zoo::rollmean, k = window, na.rm = na.rm) |>
        ts(end = end, frequency = period)
    })
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
