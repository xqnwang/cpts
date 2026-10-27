library(forecast)
library(conformalForecast)
library(dplyr)
library(tidyr)

#------------------------------------
# General setup
settings <- list(
  nobs = 5000, burnin = 500, nrep = 100, seed = 2026,
  ar_coef = c(0.8, -0.5), rho = 0.7, beta = 1, sigma = 1,
  tau = c(0.25, 1), horizon = 3, level = 90,
  fit_window = 500, cal_window = 500, calc_window = 500,
  lr = 0.1, KI = 2, Tg = 5000, delta = 0.01
)

Csat <- 2 / pi * (
  ceiling(log(settings$Tg) * settings$delta) - 1 / log(settings$Tg)
)

#------------------------------------
# Data simulation
simulate_arx <- function(tau, seed, settings) {
  set.seed(seed)
  n <- settings$nobs + settings$burnin + 2L
  u <- rnorm(n, sd = tau)
  epsilon <- rnorm(n, sd = settings$sigma)
  x <- y <- numeric(n)
  for (i in 3:n) {
    x[i] <- settings$rho * x[i - 1L] + u[i]
    y[i] <- sum(settings$ar_coef * y[c(i - 1L, i - 2L)]) +
      settings$beta * x[i] + epsilon[i]
  }
  keep <- seq.int(n - settings$nobs + 1L, n)
  list(y = ts(y[keep]), x = ts(x[keep]))
}

# Rolling ARX(2) estimates; only the future predictor inputs differ.
arx_forecast <- function(x, h, xreg, newxreg, information, rho) {
  xreg <- as.numeric(xreg)
  n <- length(x)
  train <- data.frame(
    y = as.numeric(x[3:n]),
    lag1 = as.numeric(x[2:(n - 1)]),
    lag2 = as.numeric(x[1:(n - 2)]),
    predictor = xreg[3:n]
  )
  model <- lm(y ~ lag1 + lag2 + predictor, data = train)
  future_x <- if (information == "Known") {
    as.numeric(newxreg)
  } else {
    rho^seq_len(h) * xreg[n]
  }
  history <- c(as.numeric(tail(x, 2)), rep(NA, h))
  for (j in seq_len(h)) {
    history[j + 2] <- predict(model, newdata = data.frame(
      lag1 = history[j + 1], lag2 = history[j],
      predictor = future_x[j]
    ))
  }
  list(mean = history[seq_len(h) + 2], model = model)
}

run_acmcp <- function(data, information, settings) {
  forecastfun <- function(x, h, level, xreg, newxreg) {
    arx_forecast(x, h, xreg, newxreg, information, settings$rho)
  }
  predictors <- matrix(c(data$x, rep(NA_real_, settings$horizon - 1L)),
                       ncol = 1, dimnames = list(NULL, "predictor"))
  fc <- cvforecast(
    data$y, forecastfun = forecastfun,
    h = settings$horizon, level = NULL, forward = FALSE,
    xreg = predictors, window = settings$fit_window
  )
  # Failed base fits must not be silently excluded from evaluation.
  for (h in seq_len(settings$horizon)) {
    stopifnot(all(is.finite(fc$ERROR[h:nrow(fc$ERROR), h])))
  }
  fc$level <- settings$level
  acmcp(
    fc, alpha = 1 - 0.01 * settings$level,
    ncal = settings$cal_window, rolling = TRUE,
    integrate = TRUE, scorecast = TRUE,
    lr = settings$lr, KI = settings$KI, Csat = Csat
  )
}

#------------------------------------
# Evaluation on common target dates
evaluate_acmcp <- function(object, data, settings) {
  key <- paste0(settings$level, "%")
  lower <- object$LOWER[[key]]
  upper <- object$UPPER[[key]]
  dates <- as.integer(round(time(lower)))
  first <- settings$fit_window + settings$cal_window +
    2L * settings$horizon
  rows <- which(dates >= first & dates <= settings$nobs)
  stopifnot(length(rows) >= settings$calc_window)
  target <- as.numeric(data$y)[dates[rows]]
  bind_rows(lapply(seq_len(settings$horizon), function(h) {
    lo <- as.numeric(lower[rows, h])
    up <- as.numeric(upper[rows, h])
    stopifnot(!anyNA(lo), !anyNA(up))
    covered <- target >= lo & target <= up
    # Empty intervals have zero width; unbounded intervals are retained.
    widths <- ifelse(lo > up, 0, up - lo)
    stopifnot(!anyNA(widths))
    local <- as.numeric(stats::filter(
      as.numeric(covered), rep(1 / settings$calc_window, settings$calc_window),
      sides = 1
    ))
    data.frame(
      horizon = h,
      coverage = 100 * mean(covered),
      mean_width = mean(widths),
      local_mae = 100 * mean(abs(local - settings$level / 100), na.rm = TRUE)
    )
  }))
}

#------------------------------------
# Paired repetitions
results <- list()
for (rep in seq_len(settings$nrep)) {
  for (tau in settings$tau) {
    data <- simulate_arx(tau, settings$seed + rep, settings)
    for (information in c("Known", "Forecast")) {
      message(sprintf("Replication %d/%d, tau = %g, %s",
                      rep, settings$nrep, tau, information))
      object <- run_acmcp(data, information, settings)
      results[[length(results) + 1L]] <- evaluate_acmcp(object, data, settings) |>
        mutate(replication = rep, tau = tau, information = information)
    }
  }
  saveRDS(list(
    settings = settings, Csat = Csat, results = bind_rows(results),
    completed_replications = rep, session = sessionInfo()
  ), "result/predictor_results.rds")
}
results <- bind_rows(results)

#------------------------------------
# Summary and Monte Carlo standard errors
mcse <- function(x) {
  if (length(x) < 2L || any(!is.finite(x))) return(NA_real_)
  sd(x) / sqrt(length(x))
}
metrics <- c("coverage", "mean_width", "local_mae")
# Coverage is in percent; local MAE is in percentage points.
summary <- results |>
  group_by(tau, horizon, information) |>
  summarise(
    across(all_of(metrics), list(mean = mean, mcse = mcse)),
    .groups = "drop"
  )

table <- summary |>
  mutate(
    information = factor(information, levels = c("Known", "Forecast")),
    coverage = sprintf("%.2f (%.3f)", coverage_mean, coverage_mcse),
    mean_width = sprintf("%.2f (%.3f)", mean_width_mean, mean_width_mcse),
    local_mae = sprintf("%.2f (%.3f)", local_mae_mean, local_mae_mcse)
  ) |>
  select(tau, horizon, information, all_of(metrics)) |>
  pivot_wider(names_from = information, values_from = all_of(metrics),
              names_sort = TRUE) |>
  arrange(tau, horizon)

write.csv(table, "result/predictor_table.csv",
          row.names = FALSE)
