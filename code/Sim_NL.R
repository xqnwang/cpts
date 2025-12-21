library(ggplot2)
library(forecast)
library(conformalForecast)

#------------------------------------
# General setup
## data simulation
nobs <- 2000
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
calc_window <- 100

#------------------------------------
# Data simulation from a nonlinear autoregressive process
set.seed(0)
nburnin <- 100
n <- nburnin + nobs + horizon
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
e <- rnorm(n, 0, 0.1)
y <- numeric(n)
for (t in 3:n) {
  y[t] <- sin(y[t-1]) + 0.5 * log(y[t-2] + 1) + 0.3 * x2[t] +
    0.1 * y[t-1] * x1[t]+ e[t]
}
y <- ts(y[(nburnin+1):(nburnin+nobs)], start = 1, frequency = 1)
exog_data <- data.frame(x1 = x1, x2 = x2)[-seq_len(nburnin), ] |>
  ts(start = 1, frequency = 1)
autoplot(y) +
  labs(
    title = "Nonlinear autoregressive process",
    x = "Time",
    y = ""
  ) +
  theme_bw()

# Time series cross-validation
nnar_model <- function(x, h, level, xreg, newxreg) {
  nnetar(x, p = 2, P = 0, xreg = xreg) |>
    forecast(h = h, PI = FALSE, level = level, xreg = newxreg)
}
fc <- cvforecast(y, forecastfun = nnar_model, h = horizon, level = NULL,
                 forward = forward, xreg = exog_data, window = fit_window)
# saveRDS(fc, file = "result/NL_fc.rds")

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
Tg <- 1000
delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 0.5
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
# complete.obs <- function(x) {
#   x[complete.cases(x),]
# }
# increase.cols <- function(x) {
#   ifinc <- matrix(nrow = nrow(x), ncol = ncol(x)-1)
#   for (i in 2:ncol(x)) {
#     ifinc[,i-1] <- x[,i] >= x[,i-1]
#   }
#   sum(apply(ifinc, 1, all))/nrow(ifinc)
# }
# info <- lapply(1:length(candidates), function(i) {
#   out_cov <- cov_list[[i]]
#   out_wid <- wid_list[[i]]
#   out_score <- score_list[[i]]
#   out_mean <- data.frame(
#     method = methods[i],
#     covmean = as.vector(out_cov$mean),
#     covmin = apply(out_cov$rollmean, 2, min, na.rm = TRUE),
#     covmax = apply(out_cov$rollmean, 2, max, na.rm = TRUE),
#     widmean = as.vector(out_wid$mean),
#     widmedian = as.vector(out_wid$median),
#     fanout = increase.cols(complete.obs(out_wid$width)),
#     winkler = as.vector(out_score[, paste0("Winkler_", level)]),
#     msis = as.vector(out_score[, paste0("MSIS_", level)])
#   ) |>
#     as_tibble() |>
#     rownames_to_column("horizon") |>
#     mutate(horizon = paste0("h=", horizon))
#   out_mean
# })
# NL_info <- do.call(bind_rows, info) |>
#   mutate(
#     method = factor(method, levels = methods),
#     covmean = round(covmean, 3)
#   ) |>
#   mutate(
#     covmean = covmean - 0.01*level,
#     covmin = covmin - 0.01*level,
#     covmax = covmax - 0.01*level
#   ) |>
#   arrange(horizon, method)
# NL_table <- NL_info |>
#   filter(method %in% c("MPI", "MPID", "AcMCP")) |>
#   select(horizon, covmean, widmean, winkler) |>
#   mutate(across(2:4, \(x) round(x, 3)))
# NL_table <- bind_cols(
#   filter(NL_table, horizon == "h=1"),
#   filter(NL_table, horizon == "h=2"),
#   filter(NL_table, horizon == "h=3")) |>
#   select(!starts_with("horizon")) |>
#   as.data.frame()
# NL_table <- data.frame(Methods = c("MPI", "MPID", "AcMCP"), NL_table)
# colnames(NL_table) <- c("Methods", rep(c("Coverage gap", "Mean width", "Winkler score"), 3))
# saveRDS(NL_table, file = "result/T_NL_winkler.rds")

#--------------------------
# Plots: Coverage and width
# cov_plot <- cov |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   geom_hline(yintercept = 0.01*level, linetype = "dashed", colour = "black") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Local coverage level",
#     colour = "Methods"
#   ) +
#   guides(colour = guide_legend(nrow = 1)) +
#   theme_bw()
# 
# wid_me_plot <- wid_me |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Mean interval width"
#   ) +
#   theme_bw()
# 
# wid_md_plot <- wid_md |>
#   filter(index >= fit_window + cal_window + calc_window) |>
#   filter(method != "MPI") |>
#   as_tsibble(index = index, key = c(horizon, method)) |>
#   mutate(method = factor(method, levels = methods)) |>
#   ggplot(aes(x = index, y = width, group = method, colour = method)) +
#   geom_line(size = 0.6, alpha = 0.9) +
#   scale_colour_manual(values = cols) +
#   # ggh4x::facet_grid2(cols = vars(horizon), scales = "free_y", independent = "y") +
#   facet_grid(cols = vars(horizon)) +
#   labs(
#     x = "Time",
#     y = "",
#     title = "Median interval width"
#   ) +
#   theme_bw()
# 
# P_NL_cov <- ggpubr::ggarrange(cov_plot, wid_me_plot, wid_md_plot,
#                               ncol = 1, nrow = 3,
#                               common.legend = TRUE, legend = "bottom")
# saveRDS(P_NL_cov, file = "result/P_NL_cov.rds")

#--------------------------
# Plots: Coverage and width as facets
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

P_NL_cov2 <- df |>
  ggplot(aes(x = index, y = yvar, group = method, colour = method)) +
  geom_line(size = 0.6, alpha = 0.9) +
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

saveRDS(P_NL_cov2, file = "result/P_NL_cov2.rds")

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

P_NL_box <- ggpubr::ggarrange(cov_boxplot, wid_boxplot,
                               ncol = 2, nrow = 1,
                               common.legend = TRUE, legend = "bottom")
saveRDS(P_NL_box, file = "result/P_NL_box.rds")

#--------------------------
# Plots: time plot
clip <- 200
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
P_NL_timeplot <- tp_data |>
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
# saveRDS(P_NL_timeplot, file = "result/P_NL_timeplot.rds")
