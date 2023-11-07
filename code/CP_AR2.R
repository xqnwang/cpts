library(ggplot2)
library(randomForest)
source("R/base.R")
source("R/ICP_split.R")
source("R/ICP_update.R")

# Simulation series AR(2)
set.seed(0)
series <- arima.sim(n = 500, list(ar = c(0.8897, -0.4858)), sd = sqrt(0.2))
autoplot(series) +
  ggtitle("AR(2) process") +
  xlab("Time") +
  ylab("") +
  theme_bw()
series <- as.numeric(series)

#------------------------------------------------------------
# ICP with fixed training and calibration sets
#------------------------------------------------------------
# Transfer series to data frame
p <- 2
i <- 1:400
foo <- function(k) c(series[k:(k+length(series)-p-1)])
dat <- sapply(1:(p+1), foo) |> data.frame() |> rev()
colnames(dat) <- c("y", "x1", "x2")
dat_tr <- dat[i, ]; x <- dat[i, -1]; y <- dat[i, 1]; n <- length(y)
dat_te <- dat[-i, ]; x0 <- dat[-i, -1]; y0 <- dat[-i, 1]; n0 <- length(y0)

#------
# SCP
out.SCP <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5)
mean(out.SCP$lo <= y0 & y0 <= out.SCP$up)

#------
# WCP - Weights estimated with logistic regression
zy <- as.factor(c(rep(0, n), rep(1, n0)))
zx <- rbind(x, x0) |> as.matrix()
obj.glm <- glm(zy ~ zx, family = "binomial")
prob.glm <- predict(obj.glm, type = "response")
wts.glm <- prob.glm / (1-prob.glm)

out.WCP.LR <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.glm)
mean(out.WCP.LR$lo <= y0 & y0 <= out.WCP.LR$up)

#------
# WCP - Weights estimated with random forest
obj.rf <- randomForest(zx, zy)
prob.rf <- predict(obj.rf, type = "prob")[,2]
prob.rf <- pmax(pmin(prob.rf, 0.99), 0.01)
wts.rf <- prob.rf / (1-prob.rf)

out.WCP.RF <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.rf)
mean(out.WCP.RF$lo <= y0 & y0 <= out.WCP.RF$up)

#------
# NexCP - Weights on the test set are all equal to 1
wts.exp <- c(0.99^(n+1-((1:n))), rep(1, n0))

out.NexCP <- ICP.split(dat_tr, x0, alpha = 0.1, rho = 0.5, w = wts.exp)
mean(out.NexCP$lo <= y0 & y0 <= out.NexCP$up)

# library(conformalInference)
# my.lm.funs = lm.funs()
# out <- conformal.pred.split(x, y,as.matrix(x0), w=wts.exp, split=1:200,
#                      train.fun=my.lm.funs$train,
#                      predict.fun=my.lm.funs$predict)
# mean(out$lo <= y0 & y0 <= out$up)


#------------------------------------------------------------
# ICP with sliding training and calibration sets
#------------------------------------------------------------
nfit <- 100
ncal <- 100
series0 <- series[(nfit+ncal+1):length(series)]
#------
# SCP
out_SCP <- ICP.update(series, nfit = 100, ncal = 100, alpha = 0.1)
mean(out_SCP$lo <= series0 & series0 <= out_SCP$up)

#------
# WCP - Weights estimated with logistic regression
out_WCP.LR <- ICP.update(series, nfit = 100, ncal = 100, alpha = 0.1, w = c(0,0,wts.glm))
mean(out_WCP.LR$lo <= series0 & series0 <= out_WCP.LR$up)

#------
# WCP - Weights estimated with random forest
out_WCP.RF <- ICP.update(series, nfit = 100, ncal = 100, alpha = 0.1, w = c(0,0,wts.rf))
mean(out_WCP.RF$lo <= series0 & series0 <= out_WCP.RF$up)

#------
# NexCP - Weights on the test set are all equal to 1
out_NexCP <- ICP.update(series, nfit = 100, ncal = 100, alpha = 0.1, w = wts.exp)
mean(out_NexCP$lo <= series0 & series0 <= out_NexCP$up)

#------
# ACP - Adaptively update alpha
out_ACP <- ICP.update(series, nfit = 100, ncal = 100, alpha = 0.1,
                      updateAlpha = TRUE, gamma = 0.005, updateMethod = "Simple")
mean(out_ACP$lo <= series0 & series0 <= out_ACP$up)

