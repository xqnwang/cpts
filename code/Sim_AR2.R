library(ggplot2)
library(patchwork)
library(forecast)
library(conformalForecast)

#------------------------------------
# General setup
## data simulation
nobs <- 5000
ar_coef <- c(0.8, -0.5)
sigma <- 1
!any(abs(polyroot(c(1, -ar_coef))) <= 1) # The model will be stationary if TRUE
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
calc_window <- 500

#------------------------------------
# Data simulation from AR(2)
set.seed(0)
data_AR2 <- arima.sim(n = nobs, list(ar = ar_coef), sd = sigma)
autoplot(data_AR2) +
  labs(
    title = "AR(2)",
    x = "Time",
    y = ""
  ) +
  theme_bw()

# Time series cross-validation
fc_ar2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level = level)
}
fc <- cvforecast(data_AR2, forecastfun = fc_ar2, h = horizon, level = level,
                 forward = forward, window = fit_window)

# MSCP
mscp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = NULL, kess = FALSE, quantiletype = type)

# MWCP
expweight <- function(n) 0.99^{n+1-(1:n)}
mwcp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = expweight, kess = TRUE, quantiletype = type)

# MACP
gamma <- 0.005
macp <- acp(fc, symmetric = symmetric, gamma = gamma,
            ncal = cal_window, rolling = rolling)

# MPID
Tg <- 5000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2
lr <- 0.1

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
acmcp <- acmcp(fc, ncal = cal_window, rolling = rolling,
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

#--------------------------
# Plots: Coverage and width
# cov_plot <- cov |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.8) +
#   scale_colour_manual(values = cols) +
#   geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "black") +
#   # ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Local coverage level"
#   ) +
#   theme_bw() +
#   theme(legend.position="none")
# 
# wid_me_plot <- wid_me |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.8) +
#   scale_colour_manual(values = cols) +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Mean interval width"
#   ) +
#   theme_bw() +
#   theme(legend.position="none")
# 
# wid_md_plot <- wid_md |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.8) +
#   scale_colour_manual(values = cols) +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Median interval width",
#     colour = "Methods"
#   ) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(colour = guide_legend(nrow = 1))
# 
# P_AR2_cov <- cov_plot / wid_me_plot / wid_md_plot
# saveRDS(P_AR2_cov, file = "result/P_AR2_cov.rds")

#--------------------------
# Plots: Coverage and width as panels
df <- bind_rows(
  cov |>
    filter(index >= fit_window + cal_window + calc_window) |>
    filter(method != "MPI") |>
    mutate(
      method = factor(method, levels = methods),
      yvar = coverage,
      panel = "Local coverage level"
    ),
  wid_me |>
    filter(index >= fit_window + cal_window + calc_window) |>
    filter(method != "MPI") |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Mean interval width"
    ),
  wid_md |>
    filter(index >= fit_window + cal_window + calc_window) |>
    filter(method != "MPI") |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Median interval width"
    )
) |>
  as_tsibble(index = index, key = c(horizon, method, panel))

P_AR2_cov2 <- df |>
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

saveRDS(P_AR2_cov2, file = "result/P_AR2_cov2.rds")

#--------------------------
# Boxplots: rolling coverage and width
cov_boxplot <- cov |>
  filter(method != "MPI") |>
  mutate(
    method = factor(method, levels = methods),
  ) |>
  ggplot(aes(x = method, y = coverage)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey", outliers = FALSE) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_grid(rows = vars(horizon)) +
  labs(
    x = "",
    y = "Rolling coverage",
  ) +
  theme_bw()

wid_boxplot <- wid |>
  filter(method != "MPI") |>
  mutate(
    method = factor(method, levels = methods),
  ) |>
  ggplot(aes(x = method, y = width)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey", outliers = FALSE) +
  geom_hline(
    data = filter(wid, horizon == "h=1"),
    aes(yintercept = wid |> filter(horizon == "h=1" & method == "AcMCP") |> pull(width) |> median(na.rm = TRUE)),
    linetype = "dashed", colour = "blue") +
  geom_hline(
    data = filter(wid, horizon == "h=2"),
    aes(yintercept = wid |> filter(horizon == "h=2" & method == "AcMCP") |> pull(width) |> median(na.rm = TRUE)),
    linetype = "dashed", colour = "blue") +
  geom_hline(
    data = filter(wid, horizon == "h=3"),
    aes(yintercept = wid |> filter(horizon == "h=3" & method == "AcMCP") |> pull(width) |> median(na.rm = TRUE)),
    linetype = "dashed", colour = "blue") +
  coord_flip() +
  facet_grid(rows = vars(horizon)) +
  labs(
    x = "",
    y = "Interval width",
  ) +
  theme_bw()

P_AR2_box <- ggpubr::ggarrange(cov_boxplot, wid_boxplot,
                               ncol = 2, nrow = 1,
                               common.legend = TRUE, legend = "bottom")
saveRDS(P_AR2_box, file = "result/P_AR2_box.rds")

#--------------------------
# Plots: time plot
clip <- 500
clip_start <- fit_window + cal_window + 1
clip_end <- clip_start + clip
for (h in seq_len(horizon)) {
  tb <- tibble(
    time = clip_start:clip_end,
    x = window(fc$x, start = clip_start, end = clip_end),
    horizon = paste0("h=", h),
    lower_acmcp = window(acmcp$LOWER[[1]][,h], start = clip_start, end = clip_end),
    upper_acmcp = window(acmcp$UPPER[[1]][,h], start = clip_start, end = clip_end),
    lower_mpi = window(mpi$LOWER[[1]][,h], start = clip_start, end = clip_end),
    upper_mpi = window(mpi$UPPER[[1]][,h], start = clip_start, end = clip_end),
    lower_macp = window(macp$LOWER[[1]][,h], start = clip_start, end = clip_end),
    upper_macp = window(macp$UPPER[[1]][,h], start = clip_start, end = clip_end),
  )
  assign(paste0("tp_h", h), tb)
}
tp_data <- bind_rows(mget(paste0("tp_h", seq_len(horizon))))
P_AR2_timeplot <- tp_data |>
  ggplot(aes(x = time)) +
  geom_line(aes(y = x, colour = "Observation")) +
  geom_line(aes(y = lower_macp, colour="MACP"), alpha = 1) +
  geom_line(aes(y = upper_macp, colour="MACP"), alpha = 1) +
  geom_line(aes(y = lower_mpi, colour="MPI"), alpha = 0.8) +
  geom_line(aes(y = upper_mpi, colour="MPI"), alpha = 0.8) +
  geom_line(aes(y = lower_acmcp, colour="AcMCP"), alpha = 0.7) +
  geom_line(aes(y = upper_acmcp, colour="AcMCP"), alpha = 0.7) +
  facet_grid(horizon~.) +
  scale_color_manual(
    name = "Methods",
    values = c("Observation" = "black", "MACP" = "green", "MPI" = "blue", "AcMCP" = "red"),
    breaks = c('Observation','MACP','MPI','AcMCP'),
    labels = c('Observation','MACP','MPI','AcMCP')) +
  xlab("Time") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom")
saveRDS(P_AR2_timeplot, file = "result/P_AR2_timeplot.rds")


