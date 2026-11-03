# Clip infinite MACP interval bounds.
#
# The asymmetric ACI update can push a tail's working miscoverage level to or
# below 1/(ncal + 1), at which point the conformal quantile (computed with an
# appended Inf) is infinite. Following the paper, an infinite upper (lower)
# bound for y_{t+h} is replaced by the point forecast plus (minus) the largest
# upper-tail (lower-tail) score observed up to the forecast origin t, i.e.
# max_{s <= t} e_{s|s-h} (respectively max_{s <= t} -e_{s|s-h}).
#
# In conformalForecast objects, MEAN[r, h] is the h-step forecast of y at
# target row r and ERROR[r, h] = y_r - MEAN[r, h], both indexed by target time
# from the same start. The scores available when forecasting row r at horizon
# h are therefore ERROR[1:(r - h), h].
clip_acp <- function(object) {
  pf <- as.matrix(object$MEAN)
  err <- as.matrix(object$ERROR)
  n <- nrow(pf)
  running_max <- function(x) {
    x[is.na(x)] <- -Inf
    out <- cummax(x)
    out[is.infinite(out)] <- NA_real_
    out
  }
  for (h in seq_len(ncol(pf))) {
    e <- c(err[, h], rep(NA_real_, n - nrow(err)))
    # Shift by h rows so that row r sees only scores up to row r - h
    max_up <- c(rep(NA_real_, h), running_max(e))[seq_len(n)]
    max_lo <- c(rep(NA_real_, h), running_max(-e))[seq_len(n)]
    for (k in names(object$UPPER)) {
      up <- object$UPPER[[k]]
      idx <- which(is.infinite(up[, h]))
      up[idx, h] <- pf[idx, h] + max_up[idx]
      object$UPPER[[k]] <- up
      lo <- object$LOWER[[k]]
      idx <- which(is.infinite(lo[, h]))
      lo[idx, h] <- pf[idx, h] - max_lo[idx]
      object$LOWER[[k]] <- lo
    }
  }
  object
}
