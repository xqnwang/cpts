library(ggplot2)
library(patchwork)
library(forecast)
library(conformalForecast)
reticulate::use_virtualenv("r-spci", required = TRUE)
source("code/spci.R")

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
# saveRDS(fc, file = "result/AR2_fc.rds")
# fc <- readRDS("result/AR2_fc.rds")

# MSCP
mscp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = NULL, kess = FALSE, quantiletype = type)

# MWCP
expweight <- function(n) 0.99^{n+1-(1:n)}
mwcp <- scp(fc, symmetric = symmetric, ncal = cal_window, rolling = rolling,
            weightfun = expweight, kess = TRUE, quantiletype = type)

# SPCI
spci <- spci(fc, ncal = cal_window, rolling = rolling, w = 10,
             random_state = 0L)

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

saveRDS(mget(c("mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")),
        file = "result/AR2_cp.rds")

#------------------------------------
# Result analysis
library(dplyr)
library(tidyr)
library(tibble)
library(tsibble)

candidates <- c("fc", "mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")
methods <- c("AR", "MSCP", "MWCP", "SPCI", "MACP", "MPI", "MPID", "AcMCP")

cov_list <- sapply(candidates, function(i) {
  coverage(get(i), window = calc_window, level = level)
}, simplify = FALSE,USE.NAMES = TRUE)
wid_list <- sapply(candidates, function(i) {
  width(get(i), window = calc_window, level = level, includemedian = TRUE)
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
    )
) |>
  as_tsibble(index = index, key = c(horizon, method, panel))

cols <- c(
  "AR" = "#222222", "MSCP" = "#AAAAAA", "MWCP" = "#A879CF",
  "SPCI" = "#377EB8", "MACP" = "#20AF35", "MPID" = "#F2BC12",
  "AcMCP" = "#D7191C"
)
linetypes <- c(
  "AR" = "solid", "MSCP" = "solid", "MWCP" = "dashed",
  "SPCI" = "dashed", "MACP" = "solid", "MPID" = "solid",
  "AcMCP" = "solid"
)
plot_methods <- methods[methods != "MPI"]
linewidths <- c("AR" = 0.35, "MSCP" = 0.35, "MWCP" = 0.35,
                  "SPCI" = 0.35, "MACP" = 0.45, "MPID" = 0.45, "AcMCP" = 0.45)
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

cov_plot <- df |>
  filter(panel == "Local coverage level") |>
  ggplot(aes(x = index, y = yvar, colour = method, linetype = method)) +
  geom_hline(yintercept = 0.01 * level, colour = "black",
             linetype = "dashed", linewidth = 0.35) +
  geom_line(aes(linewidth = method)) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  scale_colour_manual(values = cols, breaks = plot_methods) +
  scale_linetype_manual(values = linetypes, breaks = plot_methods) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  facet_grid(cols = vars(horizon)) +
  labs(x = NULL, y = "Rolling coverage") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  theme_cov +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

wid_me_plot <- df |>
  filter(panel == "Mean interval width") |>
  ggplot(aes(x = index, y = yvar, colour = method, linetype = method)) +
  geom_line(aes(linewidth = method)) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  scale_colour_manual(values = cols, breaks = plot_methods) +
  scale_linetype_manual(values = linetypes, breaks = plot_methods) +
  facet_grid(cols = vars(horizon)) +
  labs(x = "Time", y = "Mean interval width") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  theme_cov +
  theme(strip.text = element_blank())

P_AR2_cov2 <- (cov_plot / wid_me_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(4, 0, 0, 0),
        legend.box.spacing = grid::unit(0, "pt"))
saveRDS(P_AR2_cov2, file = "result/P_AR2_cov2.rds")

#--------------------------
# Boxplots: rolling coverage and width
cov_boxplot <- cov |>
  filter(method != "MPI") |>
  mutate(
    method = factor(method, levels = rev(methods)),
  ) |>
  ggplot(aes(x = method, y = coverage, colour = method)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey", outliers = FALSE) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "black") +
  coord_flip() +
  facet_grid(rows = vars(horizon), switch = "y") +
  labs(
    x = "",
    y = "Rolling coverage",
  ) +
  theme_cov +
  theme(legend.position = "none", strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0))

wid_boxplot <- wid |>
  filter(method != "MPI") |>
  mutate(
    method = factor(method, levels = rev(methods)),
  ) |>
  ggplot(aes(x = method, y = width, colour = method)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "grey", outliers = FALSE) +
  scale_colour_manual(values = cols) +
  geom_hline(
    data = wid |> filter(method == "AcMCP") |>
      group_by(horizon) |> summarise(width = median(width, na.rm = TRUE)),
    aes(yintercept = width),
    linetype = "dashed", colour = "black") +
  coord_flip() +
  facet_grid(rows = vars(horizon), switch = "y") +
  labs(
    x = "",
    y = "Interval width",
  ) +
  theme_cov +
  theme(legend.position = "none", strip.placement = "outside",
        strip.text.y = element_blank())

P_AR2_box <- cov_boxplot | wid_boxplot
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
  facet_grid(rows = vars(horizon), switch = "y") +
  scale_color_manual(
    name = "Methods",
    values = c("Observation" = "black", "MPI" = "blue", cols[c("MACP", "AcMCP")]),
    breaks = c('Observation','MACP','MPI','AcMCP'),
    labels = c('Observation','MACP','MPI','AcMCP')) +
  xlab("Time") +
  ylab("") +
  theme_cov +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0))
saveRDS(P_AR2_timeplot, file = "result/P_AR2_timeplot.rds")


