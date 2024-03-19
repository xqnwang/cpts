library(forecast)
source("R/lagmatrix.R")
source("R/CVforecast.R")
source("R/SCP.R")
source("R/ACP.R")
source("R/PID.R")
source("R/MCP.R")
source("R/coverage.R")
source("R/width.R")

# Simulation series AR(2)
set.seed(0)
series <- arima.sim(n = 5000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
series <- as.numeric(series)

# Setup
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |> 
    forecast(h = h, level)
}
h <- 3
level <- c(80, 95)
window <- 500
ncal <- 500
rolling <- TRUE
symmetric <- TRUE
type <- 1

fc <- CVforecast(series, forecastfun = far2, h = 3, level = level,
                 forward = TRUE, window = window, initial = 1)
coverage(fc)$mean
width(fc)$mean

## SCP
scpfc <- SCP(fc, symmetric = symmetric,
             ncal = ncal, rolling = rolling,
             weightfun = NULL, kess = FALSE,
             quantiletype = type)
coverage(scpfc)$mean
width(scpfc)$mean

## ACP
gamma <- 0.005
acpfc <- ACP(fc, symmetric = !symmetric, gamma = gamma,
             ncal = ncal, rolling = rolling,
             quantiletype = type)
coverage(acpfc)$mean
width(acpfc)$mean

## PID
Tg <- 1000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
lr <- 0.1
KI <- 2
naivefun <- function(x, h) {
  naive(x) |> 
    forecast(h = h)
}
pidfc_nsf <- PID(fc, symmetric = !symmetric,
                 ncal = ncal, rolling = rolling,
                 integrate = TRUE, scorecast = FALSE,
                 lr = lr, KI = KI, Csat = Csat)
coverage(pidfc_nsf)$mean
width(pidfc_nsf)$mean

pidfc <- PID(fc, symmetric = symmetric,
             ncal = ncal, rolling = rolling,
             integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
             lr = lr, KI = KI, Csat = Csat)
coverage(pidfc)$mean
width(pidfc)$mean

## MCP
mcpfc_nsf <- MCP(fc,
                 ncal = ncal, rolling = rolling,
                 integrate = TRUE, scorecast = FALSE,
                 lr = lr, KI = KI, Csat = Csat)
coverage(mcpfc_nsf)$mean
width(mcpfc_nsf)$mean

mcpfc <- MCP(fc,
             ncal = ncal, rolling = rolling,
             integrate = TRUE, scorecast = TRUE,
             lr = lr, KI = KI, Csat = Csat)
coverage(mcpfc)$mean
width(mcpfc)$mean

