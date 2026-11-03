library(conformalForecast)
source("code/clip_acp.R")
library(forecast)
library(tsibble)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
reticulate::use_virtualenv("./.venv", required = TRUE)
source("code/spci.R")

# Cafes, restaurants and takeaway food services in Victoria, Australia
cafe_vic <- readabs::read_abs_series("A3349417W", release_date = "2024-04-01") |>
  mutate(
    Month = yearmonth(date),
    Turnover = value
  ) |>
  filter(yearmonth(Month) <= yearmonth("2019 Dec")) |>
  as_tsibble(index = Month)
cafe <- ts(cafe_vic |> pull(Turnover), start = c(1982, 4), frequency = 12)

P_cafe_data <- ggplot(cafe_vic, aes(x = Month, y = Turnover)) +
  geom_line() +
  labs(
    x = "Year-Month",
    y = "Turnover",
    title = ""
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.25)
  )
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
type <- 1
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
saveRDS(fc, file = "result/cafe_fc.rds")
# fc <- readRDS("result/cafe_fc.rds")

# MSCP
mscp <- scp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = NULL, kess = FALSE, quantiletype = type)

# MWCP
expweight <- function(n) 0.99^{n+1-(1:n)}
mwcp <- scp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = expweight, kess = FALSE, quantiletype = type)

# SPCI
spci <- spci(fc, alpha = 1 - 0.01 * level,
             ncal = cal_window, rolling = rolling, w = 10,
             random_state = 0L)

# MACP
gamma <- 0.005
macp <- acp(fc, alpha = 1 - 0.01 * level,
            symmetric = symmetric, gamma = gamma,
            ncal = cal_window, rolling = rolling, quantiletype = type)
## imputing infinite intervals with the largest score seen so far
macp <- clip_acp(macp)

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

saveRDS(mget(c("mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")),
        file = "result/cafe_cp.rds")

#------------------------------------
# Result analysis
library(tibble)
y <- cafe

candidates <- c("mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")
methods <- c("MSCP", "MWCP", "SPCI", "MACP", "MPI", "MPID", "AcMCP")
cols <- c(
  "MSCP" = "#AAAAAA", "MWCP" = "#A879CF",
  "SPCI" = "#377EB8", "MACP" = "#20AF35", "MPID" = "#F2BC12",
  "AcMCP" = "#D7191C", "MPI" = "blue"
)
linetypes <- c(
  "MSCP" = "solid", "MWCP" = "dashed",
  "SPCI" = "dashed", "MACP" = "solid", "MPID" = "solid",
  "AcMCP" = "solid", "MPI" = "solid"
)
plot_methods <- methods
linewidths <- c("MSCP" = 0.35, "MWCP" = 0.35,
                  "SPCI" = 0.35, "MACP" = 0.45, "MPID" = 0.45, "AcMCP" = 0.45, "MPI" = 0.35)
theme_cov <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
    panel.border = element_rect(colour = "grey50", linewidth = 0.35),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text = element_text(colour = "grey20"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0, "pt"),
    panel.spacing.x = grid::unit(0.7, "lines")
  )

cov_list <- sapply(candidates, function(i) {
  conformalForecast::coverage(get(i), window = calc_window, level = level)
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
  ggplot(aes(x = horizon, y = covmean, colour = method, linetype = method)) +
  geom_line(aes(linewidth = method)) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  scale_linetype_manual(name = NULL, values = linetypes, breaks = methods) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  scale_colour_manual(name = NULL, values = cols, breaks = methods) +
  expand_limits(y = c(-0.1, 0.1)) +
  labs(
    x = "Forecast horizon",
    y = "Coverage gap",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_cov +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_width <- cafe_info |>
  mutate(horizon = str_sub(horizon, start = 3, end = -1) |> as.numeric()) |>
  ggplot(aes(x = horizon, y = widmean, colour = method, linetype = method)) +
  geom_line(aes(linewidth = method)) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  scale_linetype_manual(name = NULL, values = linetypes, breaks = methods) +
  scale_colour_manual(name = NULL, values = cols, breaks = methods) +
  labs(
    x = "Forecast horizon",
    y = "Mean interval width",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_cov +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_winkler <- cafe_info |>
  mutate(horizon = str_sub(horizon, start = 3, end = -1) |> as.numeric()) |>
  ggplot(aes(x = horizon, y = winkler, colour = method, linetype = method)) +
  geom_line(aes(linewidth = method)) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  scale_linetype_manual(name = NULL, values = linetypes, breaks = methods) +
  scale_colour_manual(name = NULL, values = cols, breaks = methods) +
  labs(
    x = "Forecast horizon",
    y = "Winkler score",
    colour = "Methods"
  ) +
  scale_x_continuous(breaks = seq_len(horizon)) +
  theme_cov +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
  # + guides(colour = guide_legend(nrow = 2, keyheight = 0.5))

P_cafe_result <- (P_cafe_covdiff / P_cafe_width / P_cafe_winkler) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(4, 0, 0, 0))
saveRDS(P_cafe_result, file = "result/P_cafe_result.rds")

P_cafe_cov <- (P_cafe_covdiff | P_cafe_width) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(0, "pt"))
saveRDS(P_cafe_cov, file = "result/P_cafe_cov.rds")
