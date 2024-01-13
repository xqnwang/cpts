library(ggplot2)
library(randomForest)
source("R/base.R")
source("R/ICP_split.R")
source("R/ICP_update.R")
source("R/PID_update.R")

# Simulation series AR(2)
set.seed(0)
series <- arima.sim(n = 5000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
autoplot(series) +
  ggtitle("AR(2) process") +
  xlab("Time") +
  ylab("") +
  theme_bw()
save(series, file = "data/AR2data.rda")
series <- as.numeric(series)

# Setup for weighted quantiles
kess <- FALSE
type <- 1

#------------------------------------------------------------
# ICP with fixed training and calibration sets
#------------------------------------------------------------
# Transfer series to data frame
p <- 2
i <- 1:1000
foo <- function(k) c(series[k:(k+length(series)-p-1)])
dat <- sapply(1:(p+1), foo) |> data.frame() |> rev()
colnames(dat) <- c("y", "x1", "x2")
dat_tr <- dat[i, ]; x <- dat[i, -1]; y <- dat[i, 1]; n <- length(y)
dat_te <- dat[-i, ]; x0 <- dat[-i, -1]; y0 <- dat[-i, 1]; n0 <- length(y0)

cov.split <- len.split <- data.frame("CP" = rep(as.double(NA), n0),
                                     "WCP.LR" = rep(as.double(NA), n0),
                                     "WCP.RF" = rep(as.double(NA), n0),
                                     "NexCP" = rep(as.double(NA), n0))
#------
# SCP
out.SCP <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5,
                     kess = kess, type = type)
cov.split$CP <- (out.SCP$lo <= y0 & y0 <= out.SCP$up)
len.split$CP <- out.SCP$up - out.SCP$lo

#------
# WCP - Weights estimated with logistic regression
zy <- as.factor(c(rep(0, n), rep(1, n0)))
zx <- rbind(x, x0) |> as.matrix()
obj.glm <- glm(zy ~ zx, family = "binomial")
prob.glm <- predict(obj.glm, type = "response")
wts.glm <- prob.glm / (1-prob.glm)

out.WCP.LR <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.glm,
                        kess = kess, type = type)
cov.split$WCP.LR <- (out.WCP.LR$lo <= y0 & y0 <= out.WCP.LR$up)
len.split$WCP.LR <- out.WCP.LR$up - out.WCP.LR$lo

#------
# WCP - Weights estimated with random forest
obj.rf <- randomForest(zx, zy)
prob.rf <- predict(obj.rf, type = "prob")[,2]
prob.rf <- pmax(pmin(prob.rf, 0.99), 0.01)
wts.rf <- prob.rf / (1-prob.rf)

out.WCP.RF <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.rf,
                        kess = kess, type = type)
cov.split$WCP.RF <- (out.WCP.RF$lo <= y0 & y0 <= out.WCP.RF$up)
len.split$WCP.RF <- out.WCP.RF$up - out.WCP.RF$lo

#------
# NexCP - Weights on the test set are all equal to 1
wts.exp <- c(0.99^(n+1-((1:n))), rep(1, n0))

out.NexCP <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.exp,
                       kess = kess, type = type)
cov.split$NexCP <- (out.NexCP$lo <= y0 & y0 <= out.NexCP$up)
len.split$NexCP <- out.NexCP$up - out.NexCP$lo

save(cov.split, len.split, file = "data/AR2_split.rda")

#------------------------------------------------------------
# ICP with rolling training and calibration sets
#------------------------------------------------------------
nfit <- 500
ncal <- 500
series0 <- series[(nfit+ncal+1):length(series)]
m0 <- length(series0)

cov.update <- len.update <- data.frame("AR" = rep(as.double(NA), m0),
                                       "CP" = rep(as.double(NA), m0),
                                       "NexCP" = rep(as.double(NA), m0),
                                       "ACP" = rep(as.double(NA), m0),
                                       "PID" = rep(as.double(NA), m0))
#------
# AR
out_AR <- AR2.update(series, nfit = nfit, ncal = ncal, alpha = 0.1)
cov.update$AR <- (out_AR$lo <= series0 & series0 <= out_AR$up)
len.update$AR <- out_AR$up - out_AR$lo

#------
# SCP
out_SCP <- ICP.update(series, nfit = nfit, ncal = ncal, alpha = 0.1,
                      kess = kess, type = type)
cov.update$CP <- (out_SCP$lo <= series0 & series0 <= out_SCP$up)
len.update$CP <- out_SCP$up - out_SCP$lo

#------
# WCP - Weights estimated with logistic regression
# out_WCP.LR <- ICP.update(series, nfit = nfit, ncal = ncal, alpha = 0.1, weight = "GLM")
# mean(out_WCP.LR$lo <= series0 & series0 <= out_WCP.LR$up)

#------
# WCP - Weights estimated with random forest
# out_WCP.RF <- ICP.update(series, nfit = nfit, ncal = ncal, alpha = 0.1, weight = "RF")
# mean(out_WCP.RF$lo <= series0 & series0 <= out_WCP.RF$up)

#------
# NexCP - Weights on the test set are all equal to 1
out_NexCP <- ICP.update(series, nfit = nfit, ncal = ncal, alpha = 0.1,
                        weight = "Exp", base = 0.99,
                        kess = kess, type = type)
cov.update$NexCP <- (out_NexCP$lo <= series0 & series0 <= out_NexCP$up)
len.update$NexCP <- out_NexCP$up - out_NexCP$lo

#------
# ACP - Adaptively update alpha
out_ACP <- ICP.update(series, nfit = nfit, ncal = ncal, alpha = 0.1,
                      updateAlpha = TRUE, gamma = 0.005, updateMethod = "Simple",
                      kess = kess, type = type)
cov.update$ACP <- (out_ACP$lo <= series0 & series0 <= out_ACP$up)
len.update$ACP <- out_ACP$up - out_ACP$lo

#------
# PID - Quantile tracking + error integration + scorecasting
out_PID <- PID.update(series, nfit = nfit, nburnin = ncal, alpha = 0.1,
                      integrate = TRUE,
                      scorecast = TRUE, ncast = ncal, # use expanding window for scorecaster if ncast = NULL
                      lr = 0.1, Csat = 1, KI = 200)
cov.update$PID <- (out_PID$lo <= series0 & series0 <= out_PID$up)
len.update$PID <- out_PID$up - out_PID$lo

save(cov.update, len.update, file = "data/AR2_update.rda")


#------------------------------------------------------------
# Coverage and width
#------------------------------------------------------------
rm(list = ls())
library(zoo)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
k <- 500

#------
# Split
load("data/AR2_split.rda")
if (NCOL(cov.split) == NCOL(len.split)) m <- NCOL(cov.split)

cov.split.mean <- apply(cov.split, 2, function(x) rollmean(x = x, k = k))
cov.split.mean <- data.frame(x = 1:nrow(cov.split.mean), cov.split.mean)
cov.split.mean.df <- pivot_longer(cov.split.mean, cols = 2:(m+1))

len.split.mean <- apply(len.split, 2, function(x) rollmean(x = x, k = k))
len.split.mean <- data.frame(x = 1:nrow(len.split.mean), len.split.mean)
len.split.mean.df <- pivot_longer(len.split.mean, cols = 2:(m+1))

len.split.median <- apply(len.split, 2, function(x) rollmedian(x = x, k = k+1))
len.split.median <- data.frame(x = 1:nrow(len.split.median), len.split.median)
len.split.median.df <- pivot_longer(len.split.median, cols = 2:(m+1))

lcl.split <- ggplot(data = cov.split.mean.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  geom_hline(yintercept=apply(cov.split, 2, mean), 
             linetype='dashed', color=c('blue', 'yellow', 'orange', 'green')) +
  scale_color_manual(
    values = c(CP="blue", WCP.LR="yellow", WCP.RF="orange", NexCP="green"),
    breaks=c('CP', 'WCP.LR', 'WCP.RF', 'NexCP')) +
  labs(title = "Local Coverage Level",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

mw.split <- ggplot(data = len.split.mean.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  scale_color_manual(
    values = c(CP="blue", WCP.LR="yellow", WCP.RF="orange", NexCP="green"),
    breaks=c('CP', 'WCP.LR', 'WCP.RF', 'NexCP')) +
  labs(title = "Mean Width",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

mdw.split <- ggplot(data = len.split.median.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  scale_color_manual(
    values = c(CP="blue", WCP.LR="yellow", WCP.RF="orange", NexCP="green"),
    breaks=c('CP', 'WCP.LR', 'WCP.RF', 'NexCP')) +
  labs(title = "Median Width",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

save(lcl.split, mw.split, mdw.split, file = "data/AR2_split_plot.rda")
ggarrange(lcl.split, mw.split, mdw.split,
          ncol = 1, nrow = 3,
          common.legend = TRUE, legend = "bottom")

#------
# Update
load("data/AR2_update.rda")
if (NCOL(cov.update) == NCOL(len.update)) m <- NCOL(cov.update)

cov.update.mean <- apply(cov.update, 2, function(x) rollmean(x = x, k = k))
cov.update.mean <- data.frame(x = 1:nrow(cov.update.mean), cov.update.mean)
cov.update.mean.df <- pivot_longer(cov.update.mean, cols = 2:(m+1))

len.update.mean <- apply(len.update, 2, function(x) rollmean(x = x, k = k))
len.update.mean <- data.frame(x = 1:nrow(len.update.mean), len.update.mean)
len.update.mean.df <- pivot_longer(len.update.mean, cols = 2:(m+1))

len.update.median <- apply(len.update, 2, function(x) rollmedian(x = x, k = k+1))
len.update.median <- data.frame(x = 1:nrow(len.update.median), len.update.median)
len.update.median.df <- pivot_longer(len.update.median, cols = 2:(m+1))

lcl.update <- ggplot(data = cov.update.mean.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  geom_hline(yintercept=apply(cov.update, 2, mean), 
             linetype='dashed', color=c('gray', 'blue', 'green', 'red', "purple")) +
  scale_color_manual(
    values = c(AR="gray", CP="blue", NexCP="green", ACP="red", PID="purple"),
    breaks=c('AR', 'CP', 'NexCP', 'ACP', 'PID')) +
  labs(title = "Local Coverage Level",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

mw.update <- ggplot(data = len.update.mean.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  scale_color_manual(
    values = c(AR="gray", CP="blue", NexCP="green", ACP="red", PID="purple"),
    breaks=c('AR', 'CP', 'NexCP', 'ACP', 'PID')) +
  labs(title = "Mean Width",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

mdw.update <- ggplot(data = len.update.median.df, aes(x = x, y = value, colour = name)) +
  geom_line() +
  scale_color_manual(
    values = c(AR="gray", CP="blue", NexCP="green", ACP="red", PID="purple"),
    breaks=c('AR', 'CP', 'NexCP', 'ACP', 'PID')) +
  labs(title = "Median Width",
       x = "Time", y = "",
       fill = "Method") +
  theme(legend.position = "bottom") +
  theme_bw()

save(lcl.update, mw.update, mdw.update, file = "data/AR2_update_plot.rda")
ggarrange(lcl.update, mw.update, mdw.update,
          ncol = 1, nrow = 3,
          common.legend = TRUE, legend = "bottom")
