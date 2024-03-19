coverage <- function(object) {
  if (any(!(c("x", "lower", "upper") %in% names(object))))
    stop("x, lower, and upper are required for coverage calculation")
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
  namelevel <- names(lower)
  
  h <- NCOL(lower[[1]])
  x <- ts(matrix(rep(object$x, h), ncol = h, byrow = FALSE),
          start = start(object$x),
          frequency = frequency(object$x))
  xend <- end(x)
  freq <- frequency(x)
  
  covmat <- lapply(namelevel, function(i) {
    lo <- window(lower[[i]], start = intvstart, end = xend)
    up <- window(upper[[i]], start = intvstart, end = xend)
    xx <- window(x, start = intvstart, end = xend)
    cov <- (lo <= xx & xx <= up)
    cov <- ts(cov, start = intvstart, end = xend, frequency = freq)
    colnames(cov) <- colnames(lo)
    return(cov)
  })
  names(covmat) <- namelevel
  
  coverage <- sapply(covmat, function(cov) {
    apply(cov, 2, mean, na.rm = TRUE)
  })
  
  out <- list(
    mean = coverage,
    ifinn = covmat
  )
  return(structure(out, class = "coverage"))
}

print.coverage <- function(x, ...) {
  x$mean
}
