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
reticulate::use_virtualenv("r-spci", required = TRUE)
source("code/spci.R")

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
  facet_grid(name ~ ., scales = "free_y", switch = "y") +
  labs(
    y = ""
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
    panel.border = element_rect(colour = "grey50", linewidth = 0.35),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 10, angle = 90),
    axis.text = element_text(colour = "grey20")
  )
P_elec_sc <- vic_elec_daily |>
  ggplot(aes(x=Temperature, y=Demand)) +
  geom_point(alpha = 0.6) +
  labs(x = "Temperature", y = "Demand") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
    panel.border = element_rect(colour = "grey50", linewidth = 0.35),
    axis.text = element_text(colour = "grey20")
  )
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

saveRDS(mget(c("mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")),
        file = "result/elec_cp.rds")

#------------------------------------
# Result analysis
library(dplyr)
library(tidyr)
library(tibble)
library(tsibble)

candidates <- c("fc", "mscp", "mwcp", "spci", "macp", "mpi", "mpid", "acmcp")
methods <- c("DR", "MSCP", "MWCP", "SPCI", "MACP", "MPI", "MPID", "AcMCP")
cols <- c(
  "DR" = "#222222", "MSCP" = "#AAAAAA", "MWCP" = "#A879CF",
  "SPCI" = "#377EB8", "MACP" = "#20AF35", "MPID" = "#F2BC12",
  "AcMCP" = "#D7191C", "MPI" = "blue"
)
linetypes <- c(
  "DR" = "solid", "MSCP" = "solid", "MWCP" = "dashed",
  "SPCI" = "dashed", "MACP" = "solid", "MPID" = "solid",
  "AcMCP" = "solid", "MPI" = "solid"
)
plot_methods <- methods
linewidths <- c("DR" = 0.35, "MSCP" = 0.35, "MWCP" = 0.35,
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
  coverage(get(i), window = calc_window, level = level)
}, simplify = FALSE,USE.NAMES = TRUE)
wid_list <- sapply(candidates, function(i) {
  width(get(i), window = calc_window, level = level, includemedian = TRUE)
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
    )
) |>
  as_tsibble(index = index, key = c(horizon, method, panel))

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
  labs(x = "Time (Year 2014)", y = "Mean interval width") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  theme_cov +
  theme(strip.text = element_blank(), axis.text.x = element_text(size = 6))

P_elec_cov2 <- (cov_plot / wid_me_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(4, 0, 0, 0),
        legend.box.spacing = grid::unit(0, "pt"),
        panel.spacing.x = grid::unit(0.3, "lines"))
saveRDS(P_elec_cov2, file = "result/P_elec_cov.rds")

#--------------------------
# Boxplots: rolling coverage and width
cov_boxplot <- cov |>
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

P_elec_box <- cov_boxplot | wid_boxplot
saveRDS(P_elec_box, file = "result/P_elec_box.rds")
