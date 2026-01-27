library(ggplot2)
library(patchwork)
library(forecast)
library(conformalForecast)

#------------------------------------
# General setup
## data simulation
nobs <- 2000
## base forecasting
horizon <- 3
level <- 90
forward <- TRUE
fit_window <- 500
## conformal prediction
symmetric <- FALSE
cal_window <- 500
rolling <- TRUE
type <- 8
## analysis
calc_window <- 100

fc <- readRDS("result/NL_fc.rds")

#------------------------------------
# Conformal prediction
## AcMCP
Tg <- 1000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 0.5
lr <- 0.1
acmcp <- acmcp(fc, alpha = 1 - 0.01 * level,
               ncal = cal_window, rolling = rolling,
               integrate = TRUE, scorecast = TRUE,
               lr = lr, KI = KI, Csat = Csat)
cov_acmcp_h1 <- coverage(acmcp, window = calc_window, level = level)$rollmean[, "h=1"]
var_cov_acmcp <- var(1 - cov_acmcp_h1, na.rm = TRUE)

## MACP
for (gamma in seq(0.005, 0.5, by = 0.015)) {
  macp_i <- acp(fc, alpha = 1 - 0.01 * level,
                symmetric = symmetric, gamma = gamma,
                ncal = cal_window, rolling = rolling)
  if (any(is.infinite(macp_i$UPPER[[paste0(level, "%")]]))) {
    upper_target <- macp_i$UPPER[[paste0(level, "%")]]
    for (j in 1:horizon) {
      prev_max <- NA_real_
      for (i in 1:NROW(upper_target)) {
        if (is.finite(upper_target[i,j])) {
          prev_max <- max(prev_max, upper_target[i, j], na.rm = TRUE)
        } else if (is.infinite(upper_target[i,j])) {
          upper_target[i,j] <- prev_max
        }
      }
    }
    macp_i$UPPER[[paste0(level, "%")]] <- upper_target
  }
  cov_i <- coverage(macp_i, window = calc_window, level = level)$rollmean[, "h=1"]
  var_cov_macp <- (1- cov_i) |> var(na.rm = TRUE)
  if (var_cov_macp <= var_cov_acmcp) {
    macp <- macp_i
    break
  }
}

#------------------------------------
# Result analysis
library(dplyr)
library(tidyr)
library(tibble)
library(tsibble)

candidates <- c("macp", "acmcp")
methods <- c("MACP", "AcMCP")
cols <- c("MACP" = "green", "AcMCP" = "red")

cov_list <- sapply(candidates, function(i) {
  coverage(get(i), window = calc_window, level = level)
}, simplify = FALSE,USE.NAMES = TRUE)
wid_list <- sapply(candidates, function(i) {
  width(get(i), window = calc_window, level = level, includemedian = TRUE)
}, simplify = FALSE,USE.NAMES = TRUE)
score_list <- sapply(candidates, function(i) {
  accuracy(get(i), byhorizon = TRUE)
}, simplify = FALSE,USE.NAMES = TRUE)

## Coverage
for (i in 1:length(candidates)) {
  out_pivot <- cov_list[[i]]$rollmean |>
    as_tsibble() |>
    mutate(horizon = key, coverage = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_cov"), out_pivot)
}
cov <- bind_rows(mget(paste0(methods, "_cov")))

## Width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$width |>
    as_tsibble() |>
    mutate(horizon = key, width = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid"), out_pivot)
}
wid <- bind_rows(mget(paste0(methods, "_wid")))

## Mean width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$rollmean |>
    as_tsibble() |>
    mutate(horizon = key, width = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid_me"), out_pivot)
}
wid_me <- bind_rows(mget(paste0(methods, "_wid_me")))

## Median width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$rollmedian |>
    as_tsibble() |>
    mutate(horizon = key, width = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid_md"), out_pivot)
}
wid_md <- bind_rows(mget(paste0(methods, "_wid_md")))

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
NL_info <- do.call(bind_rows, info) |>
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

#--------------------------
# Plots: Coverage and width as panels
df <- bind_rows(
  cov |>
    filter(index >= fit_window + cal_window + calc_window) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = coverage,
      panel = "Local coverage level"
    ),
  wid_me |>
    filter(index >= fit_window + cal_window + calc_window) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Mean interval width"
    ),
  wid_md |>
    filter(index >= fit_window + cal_window + calc_window) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Median interval width"
    )
) |>
  as_tsibble(index = index, key = c(horizon, method, panel))

P_NL_param <- df |>
  ggplot(aes(x = index, y = yvar, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(panel), cols = vars(horizon), scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    colour = "Methods"
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position="bottom")

save.image(file = "result/NL_param_workspace.RData")
# load("result/NL_param_workspace.RData")
