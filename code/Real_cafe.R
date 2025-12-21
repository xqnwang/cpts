library(conformalForecast)
library(forecast)
library(tsibbledata)
library(tsibble)
library(dplyr)
library(tidyr)
library(lubridate)
library(fable)
library(ggplot2)
library(stringr)

# Cafes, restaurants and takeaway food services in Victoria, Australia
cafe_vic <- readabs::read_abs_series("A3349417W", release_date = "2024-04-01") |>
  mutate(
    Month = yearmonth(date),
    Turnover = value
  ) |>
  filter(yearmonth(Month) <= yearmonth("2019 Dec")) |>
  as_tsibble(index = Month)
cafe <- ts(cafe_vic |> pull(Turnover), start = c(1982, 4), frequency = 12)

P_cafe_data <- cafe_vic |>
  autoplot(Turnover) +
  labs(
    x = "Year-Month",
    y = "Turnover",
    title = ""
  ) +
  theme_bw()
saveRDS(P_cafe_data, file = "result/P_cafe_data.rds")

#------------------------------------
# General setup
## base forecasting
horizon <- 12; level <- 90
forward <- FALSE
fit_window <- 20*12
## conformal prediction
symmetric <- FALSE
cal_window <- 5*12; rolling <- TRUE
type <- 8
## analysis
calc_window <- 5*12

#------------------------------------
# Time series cross-validation
comb_model <- function(x, h, level) {
  ETS <- ets(x) |> forecast(h = h)
  ARIMA <- auto.arima(x, lambda = 0, biasadj = TRUE) |> forecast(h = h)
  STL <- stlf(x, lambda = 0, h = h, biasadj = TRUE)
  
  Comb <- (ETS[["mean"]] + ARIMA[["mean"]] + STL[["mean"]])/3
  return(list(mean = Comb))
}
fc <- cvforecast(cafe, forecastfun = comb_model, h = horizon, level = NULL,
                 forward = forward, window = fit_window)
# saveRDS(fc, file = "result/cafe_fc.rds")

# MSCP
mscp <- scp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = NULL, kess = FALSE, quantiletype = type)

# MWCP
expweight <- function(n) 0.99^{n+1-(1:n)}
mwcp <- scp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = expweight, kess = TRUE, quantiletype = type)

# MACP
gamma <- 0.005
macp <- acp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, gamma = gamma,
            ncal = cal_window, rolling = rolling)

# MPID
Tg <- 200; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 100
lr <- 0.1

## PID without scorecaster
mpi <- pid(fc, alpha = 1 - 0.01 * level,
           symmetric = symmetric, ncal = cal_window, rolling = rolling,
           integrate = TRUE, scorecast = FALSE,
           lr = lr, KI = KI, Csat = Csat)

## PID with theta model as scorecaster
thetafun <- function(x, h) {
  thetaf(x, h)
}
mpid <- pid(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, ncal = cal_window, rolling = rolling,
            integrate = TRUE, scorecast = TRUE, scorecastfun = thetafun,
            lr = lr, KI = KI, Csat = Csat)

# AcMCP
acmcp <- acmcp(fc, alpha = 1 - 0.01 * level,
               ncal = cal_window, rolling = rolling,
               integrate = TRUE, scorecast = TRUE,
               lr = lr, KI = KI, Csat = Csat)

#------------------------------------
# Result analysis
library(dplyr)
library(tidyr)
library(tibble)
library(tsibble)
y <- cafe

candidates <- c("mscp", "mwcp", "macp", "mpi", "mpid", "acmcp")
methods <- c("MSCP", "MWCP", "MACP", "MPI", "MPID", "AcMCP")
cols <- c(
  "MSCP" = "yellow", "MWCP" = "orange", "MACP" = "green",
  "MPI" = "blue", "MPID" = "purple", "AcMCP" = "red"
)
linetypes <- c(
  "MSCP" = "dashed", "MWCP" = "twodash", "MACP" = "dotted",
  "MPI" = "longdash", "MPID" = "twodash", "AcMCP" = "solid"
)

cov_list <- sapply(candidates, function(i) {
  coverage(get(i), window = calc_window, level = level)
}, simplify = FALSE, USE.NAMES = TRUE)
wid_list <- sapply(candidates, function(i) {
  width(get(i), window = calc_window, level = level, includemedian = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE)
score_list <- sapply(candidates, function(i) {
  accuracy(get(i), byhorizon = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE)

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
cafe_info <- do.call(bind_rows, info) |>
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

P_cafe_covdiff <- cafe_info |>
  mutate(horizon = str_sub(horizon, start = 3, end = -1) |> as.numeric()) |>
  ggplot(aes(x = horizon, y = covmean, colour = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  scale_colour_manual(values = cols) +
  ylim(c(-0.1, 0.1)) +
  labs(
    x = "Forecast horizon",
    y = "Coverage gap",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1, keyheight = 0.5))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_width <- cafe_info |>
  mutate(horizon = str_sub(horizon, start = 3, end = -1) |> as.numeric()) |>
  ggplot(aes(x = horizon, y = widmean, colour = method)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(
    x = "Forecast horizon",
    y = "Mean interval width",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1, keyheight = 0.5))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_winkler <- cafe_info |>
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
  guides(colour = guide_legend(nrow = 1, keyheight = 0.5))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_result <- ggpubr::ggarrange(P_cafe_covdiff, P_cafe_width, P_cafe_winkler,
                                   ncol = 1, nrow = 3,
                                   common.legend = TRUE, legend = "bottom")
saveRDS(P_cafe_result, file = "result/P_cafe_result.rds")

P_cafe_cov <- ggpubr::ggarrange(P_cafe_covdiff, P_cafe_width,
                                ncol = 2, nrow = 1,
                                common.legend = TRUE, legend = "bottom")
saveRDS(P_cafe_cov, file = "result/P_cafe_cov.rds")

#--------------------------
# Plots: Coverage and width
cov_plot <- cov |>
  filter(yearmonth(index) >= yearmonth("2012 Mar") & yearmonth(index) <= yearmonth("2019 Dec")) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.9) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "black") +
  ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
  labs(
    x = "Time",
    y = "",
    title = "Local coverage level",
    colour = "Methods"
  ) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw()

wid_me_plot <- wid_me |>
  filter(yearmonth(index) >= yearmonth("2012 Mar") & yearmonth(index) <= yearmonth("2019 Dec")) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.9) +
  scale_colour_manual(values = cols) +
  ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
  labs(
    x = "Time",
    y = "",
    title = "Average interval width"
  ) +
  theme_bw()

wid_md_plot <- wid_md |>
  filter(yearmonth(index) >= yearmonth("2012 Mar") & yearmonth(index) <= yearmonth("2019 Dec")) |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.9) +
  scale_colour_manual(values = cols) +
  ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
  labs(
    x = "Time",
    y = "",
    title = "Median interval width"
  ) +
  theme_bw()

P_cafe_cov <- ggpubr::ggarrange(cov_plot, wid_me_plot, wid_md_plot,
                                ncol = 1, nrow = 3,
                                common.legend = TRUE, legend = "bottom")
# saveRDS(P_cafe_cov, file = "result/P_cafe_cov.rds")

#--------------------------
# Boxplots: rolling coverage and width
cov_boxplot <- cov |>
  mutate(
    method = factor(method, levels = methods),
  ) |>
  ggplot(aes(x = method, y = coverage)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey") +
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
    y = "Width",
  ) +
  theme_bw()
P_cafe_box <- ggpubr::ggarrange(cov_boxplot, wid_boxplot,
                                ncol = 2, nrow = 1,
                                common.legend = TRUE, legend = "bottom")
# saveRDS(P_cafe_box, file = "result/P_cafe_box.rds")

#--------------------------
# Plots: time plot
month <- cafe_vic |> distinct(Month) |> pull(Month)
clip_start <- fit_window + cal_window + 1
clip_end <- length(y)
x <- subset(fc$x, start = clip_start, end = clip_end)
for (h in seq_len(horizon)) {
  tb <- tibble(
    time = month[clip_start:clip_end],
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
P_cafe_timeplot <- tp_data |>
  ggplot(aes(x = time)) +
  geom_line(aes(y = x, colour = "Observation"), linewidth = 1) +
  geom_line(aes(y = lower_macp, colour="MACP"), alpha = 1) +
  geom_line(aes(y = upper_macp, colour="MACP"), alpha = 1) +
  #geom_line(aes(y = lower_mpi, colour="MPI"), alpha = 0.8) +
  #geom_line(aes(y = upper_mpi, colour="MPI"), alpha = 0.8) +
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
# saveRDS(P_cafe_timeplot, file = "result/P_cafe_timeplot.rds")
