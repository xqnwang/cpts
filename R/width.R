width <- function(object) {
  if (any(!(c("lower", "upper") %in% names(object))))
    stop("lower, and upper are required for coverage calculation")
  lower <- lapply(object$lower, function(lo) {
    lo_ts <- ts(as.matrix(lo),
                start = start(lo),
                frequency = frequency(lo))
    nonna <- (rowSums(is.na(lo_ts)) != NCOL(lo_ts)) |> which() |> min()
    lo_ts <- subset(lo_ts, start = nonna)
    lo_ts
  })
  upper <- lapply(object$upper, function(up) {
    up_ts <- ts(as.matrix(up),
                start = start(up),
                frequency = frequency(up))
    nonna <- (rowSums(is.na(up_ts)) != NCOL(up_ts)) |> which() |> min()
    up_ts <- subset(up_ts, start = nonna)
    up_ts
  })
  intvstart <- start(lower[[1]])
  colname <- colnames(lower[[1]])
  namelevel <- names(lower)
  freq <- frequency(lower)
  
  widmat <- lapply(namelevel, function(i) {
    wid <- upper[[i]] - lower[[i]]
    wid <- ts(wid, start = intvstart, frequency = freq)
    colnames(wid) <- colname
    return(wid)
  })
  names(widmat) <- namelevel
  
  width <- sapply(widmat, function(wid) {
    apply(wid, 2, mean, na.rm = TRUE)
  })
  
  out <- list(
    mean = width,
    width = widmat
  )
  return(structure(out, class = "width"))
}

print.width <- function(x, ...) {
  x$mean
}
