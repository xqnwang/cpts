library(ggplot2)
library(forecast)
library(conformalForecast)

#------------------------------------
# General setup
## data simulation
nobs <- 1500; ar_coef <- c(0.8, -0.6, 0.4, -0.2); sigma <- 1
!any(abs(polyroot(c(1,-ar_coef))) <= 1) # The model will be stationary if TRUE
## base forecasting
horizon <- 3; level <- 95
forward <- TRUE
fit_window <- 100
## conformal prediction
symmetric <- FALSE
cal_window <- 100; rolling <- TRUE
type <- 8
## analysis
calc_window <- 100

#------------------------------------
# Data simulation from AR(3)
set.seed(123)
data_AR4 <- arima.sim(n = nobs, list(ar = ar_coef), sd = sigma)
autoplot(data_AR4) +
  labs(
    title = "AR(4)",
    x = "Time",
    y = ""
  ) +
  theme_bw()

# Time series cross-validation
fc_ar4 <- function(x, h, level) {
  Arima(x, order = c(4, 0, 0)) |> forecast(h = h, level = level)
}
fc <- cvforecast(data_AR4, forecastfun = fc_ar4, h = horizon, level = level,
                 forward = forward, window = fit_window)

# MSCP
mscp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = NULL, kess = FALSE, quantiletype = type)

# MWCP
expweight <- function(n) 0.99^{n+1-(1:n)}
mwcp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = expweight, kess = TRUE, quantiletype = type)

# MACP
gamma <- 0.01
macp <- acp(fc, symmetric = symmetric, gamma = gamma,
            ncal = cal_window, rolling = rolling)
## replace infinite upper bounds with the overall maximum value
for (h in seq_len(horizon)) {
  inft <- is.infinite(macp$UPPER[[1]][, h])
  macp$UPPER[[1]][inft, h] <- max(macp$UPPER[[1]][!inft, h], na.rm = TRUE)
}

# MPID
Tg <- 1000; delta <- 0.001
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 3
lr <- 0.05

## PID without scorecaster
mpi <- pid(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
           integrate = TRUE, scorecast = FALSE,
           lr = lr, KI = KI, Csat = Csat)

## PID with theta model as scorecaster
thetafun <- function(x, h) {
  thetaf(x, h)
}
mpid <- pid(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            integrate = TRUE, scorecast = TRUE, scorecastfun = thetafun,
            lr = lr, KI = KI, Csat = Csat)

# AcMCP
acmcp <- mcp(fc, ncal = cal_window, rolling = rolling,
             integrate = TRUE, scorecast = TRUE,
             lr = lr, KI = KI, Csat = Csat)

#------------------------------------
# Result analysis
library(dplyr)
library(tidyr)
library(tibble)
library(tsibble)

candidates <- c("fc", "mscp", "mwcp", "macp", "mpi", "mpid", "acmcp")
methods <- c("AR", "MSCP", "MWCP", "MACP", "MPI", "MPID", "AcMCP")
cols <- c(
  "AR" = "black", "MSCP" = "yellow", "MWCP" = "orange", "MACP" = "green",
  "MPI" = "blue", "MPID" = "purple", "AcMCP" = "red"
)
linetypes <- c(
  "AR" = "solid", "MSCP" = "dashed", "MWCP" = "twodash", "MACP" = "dotted",
  "MPI" = "longdash", "MPID" = "twodash", "AcMCP" = "solid"
)

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

## Mean width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$rollmean |>
    as_tsibble() |>
    mutate(horizon = key, width = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid"), out_pivot)
}
wid <- bind_rows(mget(paste0(methods, "_wid")))

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
info <- lapply(1:length(candidates), function(i) {
  out_cov <- cov_list[[i]]
  out_wid <- wid_list[[i]]
  out_score <- score_list[[i]]
  out_mean <- data.frame(
    method = methods[i],
    covmean = as.vector(out_cov$mean),
    widmean = as.vector(out_wid$mean),
    widmedian = as.vector(out_wid$median),
    wid70 = as.vector(apply(out_wid$width, 2, quantile, probs = 0.7, na.rm = TRUE)),
    wid80 = as.vector(apply(out_wid$width, 2, quantile, probs = 0.8, na.rm = TRUE)),
    wid90 = as.vector(apply(out_wid$width, 2, quantile, probs = 0.9, na.rm = TRUE))
    # winkler = as.vector(out_score[, paste0("Winkler_", level)]),
    # msis = as.vector(out_score[, paste0("MSIS_", level)])
  ) |>
    as_tibble() |>
    rownames_to_column("horizon") |>
    mutate(horizon = paste0("h=", horizon))
  out_mean
})
info <- do.call(bind_rows, info) |>
  mutate(
    method = factor(method, levels = methods),
    covmean = round(covmean, 3)
  ) |>
  # mutate(covdiff = covmean - 0.01*level) |>
  arrange(horizon, method)

#--------------------------
# Plots: Coverage and width
cov_plot_1 <- cov |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter(!(method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "darkgray") +
  facet_grid(rows = vars(horizon)) +
  ylim(c(min(cov$coverage, na.rm = TRUE), max(cov$coverage, na.rm = TRUE))) +
  labs(
    x = "Time",
    y = "",
    title = "Local coverage level",
    colour = "Methods"
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw()
cov_plot_2 <- cov |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter((method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "darkgray") +
  facet_grid(rows = vars(horizon)) +
  ylim(c(min(cov$coverage, na.rm = TRUE), max(cov$coverage, na.rm = TRUE))) +
  labs(
    x = "Time",
    y = "",
    title = "Local coverage level",
    colour = "Methods"
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw()

wid_plot_1 <- wid |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter(!(method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(horizon), scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title = "Average interval width"
  ) +
  theme_bw()
wid_plot_2 <- wid |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter((method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(horizon), scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title = "Average interval width"
  ) +
  theme_bw()

wid_md_plot_1 <- wid_md |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter(!(method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(horizon), scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title = "Median interval width"
  ) +
  theme_bw()
wid_md_plot_2 <- wid_md |>
  filter(index >= fit_window + cal_window + calc_window) |>
  filter((method %in% c("MPI","MPID","AcMCP"))) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(horizon), scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title = "Median interval width"
  ) +
  theme_bw()

p_part1 <- ggpubr::ggarrange(cov_plot_1, wid_plot_1, wid_md_plot_1,
                             ncol = 3, nrow = 1,
                             common.legend = TRUE, legend = "bottom")
p_part2 <- ggpubr::ggarrange(cov_plot_2, wid_plot_2, wid_md_plot_2,
                             ncol = 3, nrow = 1,
                             common.legend = TRUE, legend = "bottom")
ggpubr::ggarrange(p_part1, p_part2, nrow = 2)

#--------------------------
# Plots: quantile estimate of forecast errors
corplot_error <- function(object, xindex, yindex, label, upper) {
  x <- object[complete.cases(object), c(xindex, yindex)]
  colnames(x) <- c("V1", "V2")
  ggplot(as.data.frame(x), aes(x = V1, y = V2)) +
    geom_point() +
    labs(
      x = paste0("h=", xindex),
      y = paste0("h=", yindex),
      title = paste0(label, " (", sprintf("%.2f", cor(x[,1], x[,2])), ") - ",
                     ifelse(upper, "Upper", "Lower"))
    ) +
    theme_bw()
}
corplot <- function(object, xindex, yindex, label, upper = TRUE) {
  if (upper) {
    x <- (object$UPPER[[1]]-object$MEAN)[complete.cases(object$UPPER[[1]]-object$MEAN),c(xindex, yindex)]
  } else {
    x <- (object$LOWER[[1]]-object$MEAN)[complete.cases(object$LOWER[[1]]-object$MEAN),c(xindex, yindex)]
  }
  colnames(x) <- c("V1", "V2")
  ggplot(as.data.frame(x), aes(x = V1, y = V2)) +
    geom_point() +
    labs(
      x = paste0("h=", xindex),
      y = paste0("h=", yindex),
      title = paste0(label, " (", sprintf("%.2f", cor(x[,1], x[,2])), ")")
    ) +
    theme_bw()
}

library(patchwork)
for (ifupper in c(TRUE, FALSE)) {
  p1 <- corplot_error(fc$ERROR, 1, 2, "Error", ifupper)
  p2 <- corplot_error(fc$ERROR, 1, 3, "Error", ifupper)
  p3 <- corplot_error(fc$ERROR, 2, 3, "Error", ifupper)
  
  p10 <- corplot(fc, 1, 2, "AR", ifupper)
  p20 <- corplot(fc, 1, 3, "AR", ifupper)
  p30 <- corplot(fc, 2, 3, "AR", ifupper)
  
  p11 <- corplot(mscp, 1, 2, "MSCP", ifupper)
  p21 <- corplot(mscp, 1, 3, "MSCP", ifupper)
  p31 <- corplot(mscp, 2, 3, "MSCP", ifupper)
  
  p12 <- corplot(mscp, 1, 2, "MWCP", ifupper)
  p22 <- corplot(mscp, 1, 3, "MWCP", ifupper)
  p32 <- corplot(mscp, 2, 3, "MWCP", ifupper)
  
  p13 <- corplot(macp, 1, 2, "MACP", ifupper)
  p23 <- corplot(macp, 1, 3, "MACP", ifupper)
  p33 <- corplot(macp, 2, 3, "MACP", ifupper)
  
  p14 <- corplot(mpi, 1, 2, "MPI", ifupper)
  p24 <- corplot(mpi, 1, 3, "MPI", ifupper)
  p34 <- corplot(mpi, 2, 3, "MPI", ifupper)
  
  p15 <- corplot(mpid, 1, 2, "MPID", ifupper)
  p25 <- corplot(mpid, 1, 3, "MPID", ifupper)
  p35 <- corplot(mpid, 2, 3, "MPID", ifupper)
  
  p16 <- corplot(acmcp, 1, 2, "AcMCP", ifupper)
  p26 <- corplot(acmcp, 1, 3, "AcMCP", ifupper)
  p36 <- corplot(acmcp, 2, 3, "AcMCP", ifupper)
  pp <- (p1 | p2 | p3) / (p10 | p20 | p30) / (p11 | p21 | p31) / (p12 | p22 | p32) /  (p13 | p23 | p33) / (p14 | p24 | p34) / (p15 | p25 | p35) / (p16 | p26 | p36)
  print(pp)
}

#--------------------------
# Plots: forecast errors
lo <- `colnames<-` ((acmcp$LOWER[[1]] - acmcp$MEAN),
                    paste0("h=",seq_len(horizon))) |>
  as_tsibble() |>
  mutate( horizon = key, lower = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  as_tibble()

up <- `colnames<-` ((acmcp$UPPER[[1]] - acmcp$MEAN),
                    paste0("h=",seq_len(horizon))) |>
  as_tsibble() |>
  mutate(horizon = key, upper = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  as_tibble()

lo_mpi <- `colnames<-` ((mpi$LOWER[[1]] - mpi$MEAN),
                         paste0("h=",seq_len(horizon))) |>
  as_tsibble() |>
  mutate( horizon = key, lower_mpi = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  as_tibble()

up_mpi <- `colnames<-` ((mpi$UPPER[[1]] - mpi$MEAN),
                         paste0("h=",seq_len(horizon))) |>
  as_tsibble() |>
  mutate(horizon = key, upper_mpi = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  as_tibble()

er <- fc$ERROR |>
  as_tsibble() |>
  mutate(horizon = key, error = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  as_tibble()

full_join(lo, up) |>
  full_join(lo_mpi) |>
  full_join(up_mpi) |>
  full_join(er) |>
  pivot_longer(
    cols = lower:error,
    names_to = "line",
    values_to = "forecast"
  ) |>
  as_tsibble(index = index, key = c(horizon, line)) |>
  filter(index > fit_window+cal_window) |>
  #filter(between(index, 201, 700)) |>
  ggplot(aes(x = index, y = forecast, group = line, colour = line)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(
    values = c("error" = "black",
               "lower" = "red", "upper" = "red",
               "lower_mpi" = "blue", "upper_mpi" = "blue")
  ) +
  facet_grid(horizon~.) +
  xlab("Time") +
  ylab("Forecast error") +
  theme_bw()

#-----------------
## Others
width_all <- lapply(1:length(candidates), function(i) {
  cov <- cov_list[[i]]$mean |>
    enframe(name = "horizon", value = "coverage")
  wid_list[[i]]$width |>
    as_tibble() |>
    pivot_longer(
      everything(),
      names_to = "horizon",
      values_to = "width"
    ) |>
    mutate(method = methods[i]) |>
    filter(!is.na(width)) |>
    left_join(cov, by = "horizon")
})
width_all <- bind_rows(width_all) |>
  mutate(method = factor(method, levels = methods))

ggplot(width_all, aes(x = coverage, y = width, group = method, colour = method)) +
  geom_boxplot(width = 0.0005, position = position_dodge(width = 0),
               fill = NA, alpha = 0.5, size = 0.5, varwidth = TRUE) +
  geom_vline(xintercept = 0.01*level, linetype = "dashed", colour = "darkgray") +
  scale_x_continuous() +
  scale_colour_manual(values = cols) +
  facet_grid(cols = vars(horizon), scales = "free_x") +
  labs(
    x = "Coverage",
    y = "Interval width"
  ) +
  theme_bw()








