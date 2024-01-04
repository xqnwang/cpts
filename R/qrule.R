## Discrete quantiles (types 1-3) & Continuous quantiles (types 4-9)
weighted.quantile <- function(x, w = NULL, p,
                              qrule = c("math", "school", "shahvaish",
                                        "hf1", "hf2", "hf3",
                                        "hf4", "hf5", "hf6", "hf7", "hf8", "hf9")){
  if (is.null(w)) w <- rep(1, length(x))
  
  # Remove zero weights and corresponding elements
  if (any(zero <- (w == 0))){
    w <- w[!zero]
    x <- x[!zero]
  }
  
  if (qrule %in% c("math", "hf1")){
    qdata <- qs(x, w, p)
    if (qdata$wlow == 0) qdata$qlow else qdata$qup
  } else if (qrule %in% c("school", "hf2")){
    qdata <- qs(x, w, p)
    if (qdata$wlow == 0) (qdata$qlow + qdata$qup)/2 else qdata$qup
  } else if (qrule == "hf3"){
    # Only consider distinct data values
    w <- rowsum(w, x)
    x <- sort(unique(x))
    qdata <- qs(x,w,p)
    if ((qdata$wlow == 0) && (qdata$ilow %% 2 == 0)) qdata$qlow else qdata$qup
  } else if (qrule == "hf4"){
    qdata <- qs(x,w,p)
    gamma <- with(qdata, wlow/(wup + wlow))
    qdata$qlow*(1 - gamma) + qdata$qup*gamma
  } else if (qrule %in% c("hf5", "hf6", "hf7", "hf8", "hf9")){
    n <- length(x)
    if (n == 1) return(x)
    
    ii <- order(x)
    x <- x[ii]
    w <- w[ii]
    cumw <- cumsum(w)
    
    if (qrule == "hf5"){
      pk <- (cumw - w/2)/(cumw[n])
    } else if (qrule == "hf6"){
      pk <- cumw/(cumw[n] + w[n])
    } else if (qrule == "hf7"){
      pk <- c(0, cumw[-n])/cumw[n-1]
    } else if (qrule == "hf8"){
      pk <- (cumw - w/3)/(cumw[n] + w[n]/3)
    } else if (qrule == "hf9"){
      pk <- c(cumw - w*3/8)/(cumw[n] + w[n]/4)
    }
    approx(pk, x, p, method = "linear", rule = 2)$y
  } else if (qrule == "shahvaish"){
    n <- length(x)
    if (n == 1) return(x)
    
    ii <- order(x)
    x <- x[ii]
    w <- w[ii]
    wbar <- w/mean(w)
    S <- cumsum(wbar)
    pk <- (S + 1/2 - w/2)/(n + 1)
    approx(pk, x, p, method = "constant", f = 0, rule = 2)$y
  }
}

last <- function(x) {
  if (any(x)) max(which(x)) else 1
}

qs <- function(x, w, p){
  n <- length(x)
  ii <- order(x)
  x <- x[ii]
  cumw <- cumsum(w[ii])
  
  pos <- last(cumw <= p*sum(w))
  posnext <- if (pos == length(x)) pos else pos+1
  
  list(qlow = x[pos], qup = x[posnext],
       ilow = pos, iup = posnext,
       wlow = p - cumw[pos]/sum(w), wup = cumw[posnext]/sum(w) - p)
}
