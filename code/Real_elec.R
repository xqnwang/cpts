library(conformalForecast)
library(forecast)
library(tsibbledata)
library(tsibble)
library(dplyr)
library(tidyr)
library(lubridate)
library(fable)
library(ggplot2)
library(patchwork)
library(stringr)

vic_elec_daily <- vic_elec |>
  index_by(Date = date(Time)) |>
  summarise(
    Demand = sum(Demand) / 1e3,
    Temperature = max(Temperature),
    Holiday = any(Holiday)
  ) |>
  mutate(
    DOW = wday(Date, label = TRUE),
    Workday = !Holiday & !(DOW %in% c("Sat", "Sun")),
    Cooling = pmax(Temperature, 18)
  ) |>
  select(-c(DOW, Holiday))
P_elec_ts <- vic_elec_daily |>
  pivot_longer(c(Demand, Temperature)) |>
  ggplot(aes(x = Date, y = value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free_y") +
  labs(
    y = ""
  ) +
  theme_bw()
P_elec_sc <- vic_elec_daily |>
  ggplot(aes(x=Temperature, y=Demand)) +
  geom_point(alpha = 0.6) +
  labs(x = "Temperature", y = "Demand") +
  theme_bw()
P_elec_data <- ggpubr::ggarrange(P_elec_ts, P_elec_sc,
                                 ncol = 2, nrow = 1)
saveRDS(P_elec_data, file = "result/P_elec_data.rds")

#------------------------------------
# General setup
## base forecasting
horizon <- 7
level <- 90
forward <- TRUE
fit_window <- vic_elec_daily |> filter(year(Date) < 2014) |> nrow()
## data
elec <- ts(vic_elec_daily, start = c(1,1), frequency = 7)[,-1]
y <- elec[,1] |> head(-horizon)
exog_data <- elec[,-1]
## conformal prediction
symmetric <- FALSE
cal_window <- 100
rolling <- TRUE
type <- 8
## analysis
calc_window <- 100

#------------------------------------
# Time series cross-validation
dyreg_model <- function(x, h, level, xreg, newxreg) {
  auto.arima(x, xreg = xreg) |>
    forecast(h = h, level = level, xreg = newxreg)
}
fc <- cvforecast(y, forecastfun = dyreg_model, h = horizon, level = level,
                 forward = forward, xreg = exog_data, window = fit_window)
# saveRDS(fc, file = "result/elec_fc.rds")
# fc <- readRDS("result/elec_fc.rds")

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
## imputing infinite intervals with the largest score seen so far
if (any(is.infinite(macp$UPPER[[paste0(level, "%")]]))) {
  upper_target <- macp$UPPER[[paste0(level, "%")]]
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
  macp$UPPER[[paste0(level, "%")]] <- upper_target
}

# MPID
Tg <- 1000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 30
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
methods <- c("DR", "MSCP", "MWCP", "MACP", "MPI", "MPID", "AcMCP")
cols <- c(
  "DR" = "black", "MSCP" = "yellow", "MWCP" = "orange", "MACP" = "green",
  "MPI" = "blue", "MPID" = "purple", "AcMCP" = "red"
)
linetypes <- c(
  "DR" = "solid", "MSCP" = "dashed", "MWCP" = "twodash", "MACP" = "dotted",
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
    as_tibble() |>
    mutate(index = vic_elec_daily$Date[length(y)-n():1+1]) |>
    as_tsibble(index = index) |>
    pivot_longer(seq_len(horizon), names_to = "horizon", values_to = "coverage") |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_cov"), out_pivot)
}
cov <- bind_rows(mget(paste0(methods, "_cov")))

## Width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$width |>
    as_tibble() |>
    mutate(index = vic_elec_daily$Date[length(y)-n():1+1]) |>
    as_tsibble(index = index) |>
    pivot_longer(seq_len(horizon), names_to = "horizon", values_to = "width") |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid"), out_pivot)
}
wid <- bind_rows(mget(paste0(methods, "_wid")))

## Mean width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$rollmean |>
    as_tibble() |>
    mutate(index = vic_elec_daily$Date[length(y)-n():1+1]) |>
    as_tsibble(index = index) |>
    pivot_longer(seq_len(horizon), names_to = "horizon", values_to = "width") |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid_me"), out_pivot)
}
wid_me <- bind_rows(mget(paste0(methods, "_wid_me")))

## Median width
for (i in 1:length(candidates)) {
  out_pivot <- wid_list[[i]]$rollmedian |>
    as_tibble() |>
    mutate(index = vic_elec_daily$Date[length(y)-n():1+1]) |>
    as_tsibble(index = index) |>
    pivot_longer(seq_len(horizon), names_to = "horizon", values_to = "width") |>
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
elec_info <- do.call(bind_rows, info) |>
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

P_elec_winkler <- elec_info |>
  mutate(horizon = str_sub(horizon, start = 3, end = -1) |> as.numeric()) |>
  ggplot(aes(x = horizon, y = winkler, colour = method)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(
    x = "Forecast horizon",
    y = "Winkler score",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 2, keyheight = 0.5))
saveRDS(P_elec_winkler, file = "result/P_elec_winkler.rds")

#--------------------------
# Plots: Coverage and width
# cov_plot <- cov |>
#   filter(index >= ymd("2014-07-19")) |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "black") +
#   # ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time (Year 2014)",
#     y = "",
#     title = "Local coverage level"
#   ) +
#   theme_bw() +
#   theme(legend.position="none")
# 
# wid_me_plot <- wid_me |>
#   filter(index >= ymd("2014-07-19")) |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   # ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time (Year 2014)",
#     y = "",
#     title = "Mean interval width"
#   ) +
#   theme_bw() +
#   theme(legend.position="none")
# 
# wid_md_plot <- wid_md |>
#   filter(index >= ymd("2014-07-19")) |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   # ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time (Year 2014)",
#     y = "",
#     title = "Median interval width",
#     colour = "Methods"
#   ) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(colour = guide_legend(nrow = 1))
# 
# P_elec_cov <- cov_plot / wid_me_plot / wid_md_plot
# saveRDS(P_elec_cov, file = "result/P_elec_cov.rds")

#--------------------------
# Plots: Coverage and width as panels
df <- bind_rows(
  cov |>
    filter(index >= ymd("2014-07-19")) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = coverage,
      panel = "Local coverage level"
    ),
  wid_me |>
    filter(index >= ymd("2014-07-19")) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Mean interval width"
    ),
  wid_md |>
    filter(index >= ymd("2014-07-19")) |>
    mutate(
      method = factor(method, levels = methods),
      yvar = width,
      panel = "Median interval width"
    )
) |>
  as_tsibble(index = index, key = c(horizon, method, panel))

P_elec_cov2 <- df |>
  ggplot(aes(x = index, y = yvar, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.9) +
  scale_colour_manual(values = cols) +
  facet_grid(rows = vars(panel),cols = vars(horizon), scales = "free_y") +
  labs(
    x = "Time (Year 2014)",
    y = "",
  ) +
  theme_bw() +
  theme(legend.position="bottom")

saveRDS(P_elec_cov2, file = "result/P_elec_cov2.rds")

#--------------------------
# Boxplots: rolling coverage and width
cov_boxplot <- cov |>
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
  mutate(
    method = factor(method, levels = methods),
  ) |>
  ggplot(aes(x = method, y = width)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey", outliers = FALSE)
for (h in seq_len(horizon)) {
  wid_boxplot <- wid_boxplot +
    geom_hline(
      data = filter(wid, horizon == paste0("h=", h)),
      aes_string(yintercept = wid |> filter(horizon == paste0("h=", h) & method == "AcMCP") |> pull(width) |> median(na.rm = TRUE)),
      linetype = "dashed", colour = "blue")
}
wid_boxplot <- wid_boxplot +
  coord_flip() +
  facet_grid(rows = vars(horizon)) +
  labs(
    x = "",
    y = "Interval width",
  ) +
  theme_bw()
P_elec_box <- ggpubr::ggarrange(cov_boxplot, wid_boxplot,
                               ncol = 2, nrow = 1,
                               common.legend = TRUE, legend = "bottom")
saveRDS(P_elec_box, file = "result/P_elec_box.rds")

#--------------------------
# Plots: time plot
clip_start <- fit_window + cal_window + 1
clip_end <- length(y)
x <- subset(fc$x, start = clip_start, end = clip_end)
for (h in seq_len(horizon)) {
  tb <- tibble(
    time = (distinct(cov, index) |> pull(index))[-1],
    x = x,
    horizon = paste0("h=", h),
    lower_acmcp = window(acmcp$LOWER[[1]][,h], start = start(x), end = end(x)),
    upper_acmcp = window(acmcp$UPPER[[1]][,h], start = start(x), end = end(x)),
    lower_mpi = window(mpi$LOWER[[1]][,h], start = start(x), end = end(x)),
    upper_mpi = window(mpi$UPPER[[1]][,h], start = start(x), end = end(x)),
    lower_macp = window(macp$LOWER[[1]][,h], start = start(x), end = end(x)),
    upper_macp = window(macp$UPPER[[1]][,h], start = start(x), end = end(x)),
  )
  assign(paste0("tp_h", h), tb)
}
tp_data <- bind_rows(mget(paste0("tp_h", seq_len(horizon))))
P_elec_timeplot <- tp_data |>
  ggplot(aes(x = time)) +
  geom_line(aes(y = x, colour = "Observation"), linewidth = 1) +
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
  xlab("Time (Year 2014)") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom")
saveRDS(P_elec_timeplot, file = "result/P_elec_timeplot.rds")
