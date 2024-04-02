library(forecast)
source("R/lagmatrix.R")
source("R/cvforecast.R")
source("R/errors.R")
source("R/SCP.R")
source("R/ACP.R")
source("R/PID.R")
source("R/MCP.R")
source("R/coverage.R")
source("R/width.R")

## Simulation series AR(2)
set.seed(0)
series <- arima.sim(n = 5000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
series <- as.numeric(series)

## Setup
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level)
}
h <- 3
level <- 90
window <- 500
ncal <- 500
rolling <- TRUE
symmetric <- FALSE
type <- 1

## CV
fc <- cvforecast(series, forecastfun = far2, h = 3, level = level,
                 forward = TRUE, window = window, initial = 1)
fc_score <- accuracy(fc, byhorizon = TRUE)
fc_cov <- coverage(fc, window = 500, level = level)
fc_wid <- width(fc, window = 500, level = level, includemedian = TRUE)

## SCP
scpfc <- scp(fc, symmetric = symmetric,
             ncal = ncal, rolling = rolling,
             weightfun = NULL, kess = FALSE,
             quantiletype = type)
scpfc_score <- accuracy(scpfc, byhorizon = TRUE)
scpfc_cov <- coverage(scpfc, window = 500, level = level)
scpfc_wid <- width(scpfc, window = 500, level = level, includemedian = TRUE)

## ACP
gamma <- 0.005
acpfc <- acp(fc, symmetric = symmetric, gamma = gamma,
             ncal = ncal, rolling = rolling,
             quantiletype = type)
acpfc_score <- accuracy(acpfc, byhorizon = TRUE)
acpfc_cov <- coverage(acpfc, window = 500, level = level)
acpfc_wid <- width(acpfc, window = 500, level = level, includemedian = TRUE)

## PID
Tg <- 1000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
lr <- 0.05
KI <- 2
naivefun <- function(x, h) {
  naive(x) |> forecast(h = h)
}
pidfc_nsf <- pid(fc, symmetric = symmetric,
                 ncal = ncal, rolling = rolling,
                 integrate = TRUE, scorecast = FALSE,
                 lr = lr, KI = KI, Csat = Csat)
pidfc_nsf_score <- accuracy(pidfc_nsf, byhorizon = TRUE)
pidfc_nsf_cov <- coverage(pidfc_nsf, window = 500, level = level)
pidfc_nsf_wid <- width(pidfc_nsf, window = 500, level = level, includemedian = TRUE)

pidfc <- pid(fc, symmetric = symmetric,
             ncal = ncal, rolling = rolling,
             integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
             lr = lr, KI = KI, Csat = Csat)
pidfc_score <- accuracy(pidfc, byhorizon = TRUE)
pidfc_cov <- coverage(pidfc, window = 500, level = level)
pidfc_wid <- width(pidfc, window = 500, level = level, includemedian = TRUE)

## MCP
mcpfc <- mcp(fc,
             ncal = ncal, rolling = rolling,
             integrate = TRUE, scorecast = TRUE,
             lr = lr, KI = KI, Csat = Csat)
mcpfc_score <- accuracy(mcpfc, byhorizon = TRUE)
mcpfc_cov <- coverage(mcpfc, window = 500, level = level)
mcpfc_wid <- width(mcpfc, window = 500, level = level, includemedian = TRUE)

## Interval coverage and width plot
library(tsibble)
library(tidyverse)

dat_cov <- fc_cov$rolling
dat_cov |>
  as_tsibble() |>
  mutate(
    horizon = key,
    coverage = value) |>
  select(-c(key, value)) |>
  update_tsibble(key = horizon) |>
  autoplot(coverage) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed") +
  facet_grid(vars(horizon)) +
  theme(legend.position = "none")

dat_widmean <- fc_wid$rollmean
dat_widmean |>
  as_tsibble() |>
  mutate(
    horizon = key,
    mean_width = value) |>
  select(-c(key, value)) |>
  update_tsibble(key = horizon) |>
  autoplot(mean_width) +
  facet_grid(vars(horizon)) +
  theme(legend.position = "none")

dat_widmedian <- fc_wid$rollmedian
dat_widmedian |>
  as_tsibble() |>
  mutate(
    horizon = key,
    median_width = value) |>
  select(-c(key, value)) |>
  update_tsibble(key = horizon) |>
  autoplot(median_width) +
  facet_grid(vars(horizon)) +
  theme(legend.position = "none")
