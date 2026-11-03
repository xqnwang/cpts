library(forecast)
library(conformalForecast)
source("code/clip_acp.R")

#------------------------------------
# General setup
horizon <- 3
level <- 90
## conformal prediction
symmetric <- FALSE
cal_window <- 500
rolling <- TRUE
type <- 1
## analysis
calc_window <- 500

fc <- readRDS("result/AR2_fc.rds")

#------------------------------------
# Conformal prediction
## AcMCP
Tg <- 5000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2
lr <- 0.1
acmcp <- acmcp(fc, ncal = cal_window, rolling = rolling,
               integrate = TRUE, scorecast = TRUE,
               lr = lr, KI = KI, Csat = Csat)
cov_acmcp_h1 <- conformalForecast::coverage(acmcp, window = calc_window, level = level)$rollmean[, "h=1"]
var_cov_acmcp <- var(1 - cov_acmcp_h1, na.rm = TRUE)

## MACP
macp <- NULL
for (gamma in seq(0.005, 0.5, by = 0.015)) {
  macp_i <- acp(fc, symmetric = symmetric, gamma = gamma,
                ncal = cal_window, rolling = rolling, quantiletype = type)
  macp_i <- clip_acp(macp_i)
  cov_i <- conformalForecast::coverage(macp_i, window = calc_window, level = level)$rollmean[, "h=1"]
  var_cov_macp <- (1- cov_i) |> var(na.rm = TRUE)
  if (var_cov_macp <= var_cov_acmcp) {
    macp <- macp_i
    break
  }
}

if (is.null(macp)) stop("No gamma satisfies the post-hoc matching criterion.")

#------------------------------------
# Result analysis
library(dplyr)
library(tibble)

candidates <- c("macp", "acmcp")
methods <- c("MACP", "AcMCP")

cov_list <- sapply(candidates, function(i) {
  conformalForecast::coverage(get(i), window = calc_window, level = level)
}, simplify = FALSE,USE.NAMES = TRUE)
wid_list <- sapply(candidates, function(i) {
  width(get(i), window = calc_window, level = level, includemedian = TRUE)
}, simplify = FALSE,USE.NAMES = TRUE)
score_list <- sapply(candidates, function(i) {
  accuracy(get(i), byhorizon = TRUE)
}, simplify = FALSE,USE.NAMES = TRUE)

## Overall measurements
complete.obs <- function(x) {
  x[complete.cases(x),]
}
increase.cols <- function(x) {
  ifinc <- matrix(nrow = nrow(x), ncol = ncol(x)-1)
  for (i in 2:ncol(x)) {
    ifinc[,i-1] <- x[,i] >= x[,i-1]
  }
  sum(apply(ifinc, 1, all))/nrow(ifinc)
}
info <- lapply(1:length(candidates), function(i) {
  out_cov <- cov_list[[i]]
  out_wid <- wid_list[[i]]
  out_score <- score_list[[i]]
  out_mean <- data.frame(
    method = methods[i],
    covmean = as.vector(out_cov$mean),
    covmin = apply(out_cov$rollmean, 2, min, na.rm = TRUE),
    covmax = apply(out_cov$rollmean, 2, max, na.rm = TRUE),
    widmean = as.vector(out_wid$mean),
    widmedian = as.vector(out_wid$median),
    fanout = increase.cols(complete.obs(out_wid$width)),
    winkler = as.vector(out_score[, paste0("Winkler_", level)]),
    msis = as.vector(out_score[, paste0("MSIS_", level)])
  ) |>
    as_tibble() |>
    rownames_to_column("horizon") |>
    mutate(horizon = paste0("h=", horizon))
  out_mean
})
AR2_info <- do.call(bind_rows, info) |>
  mutate(
    method = factor(method, levels = methods),
    covmean = round(covmean, 3)
  ) |>
  mutate(
    covmean = covmean - 0.01*level,
    covmin = covmin - 0.01*level,
    covmax = covmax - 0.01*level
  ) |>
  arrange(horizon, method)
