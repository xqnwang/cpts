# Some errors in the function, see docs/quantile.qmd for more details

## Discrete quantiles (types 1-3), Continuous quantiles (types 4-9), and Shah & Vaish
weighted.quantile.qrule <- function(x, w = NULL, p, type = 7){
  w <- w %||% rep(1, length(x))
  w <- w/sum(w)
  
  # Remove zero weights and corresponding elements
  if (any(zero <- (w == 0))){
    w <- w[!zero]
    x <- x[!zero]
  }
  
  if (type == 1){
    # the method commonly used in conformal prediction papers
    qdata <- qs(x, w, p)
    if (qdata$wlow == 0) qdata$qlow else qdata$qup
  } else if (type == 2){
    qdata <- qs(x, w, p)
    if (qdata$wlow == 0) (qdata$qlow + qdata$qup)/2 else qdata$qup
  } else if (type == 3){
    # Only consider distinct data values
    w <- rowsum(w, x)
    x <- sort(unique(x))
    qdata <- qs(x, w, p)
    if ((qdata$wlow == 0) && (qdata$ilow %% 2 == 0)) qdata$qlow else qdata$qup
  } else if (type == 4){
    qdata <- qs(x,w,p)
    gamma <- with(qdata, wlow/(wup + wlow))
    qdata$qlow*(1 - gamma) + qdata$qup*gamma
  } else if (type %in% 5:9){
    n <- length(x)
    if (n == 1) return(x)
    
    ii <- order(x)
    x <- x[ii]
    w <- w[ii]
    cumw <- cumsum(w)
    
    pk = switch(type - 4,
                # type 5
                (cumw - w/2)/(cumw[n]),
                # type 6
                cumw/(cumw[n] + w[n]),
                # type 7
                c(0, cumw[-n])/cumw[n-1],
                # type 8
                (cumw - w/3)/(cumw[n] + w[n]/3),
                # type 9
                c(cumw - w*3/8)/(cumw[n] + w[n]/4)
    )
    approx(pk, x, p, method = "linear", rule = 2)$y
  } 
  # else if (type == "shahvaish"){
  #   n <- length(x)
  #   if (n == 1) return(x)
  #   
  #   ii <- order(x)
  #   x <- x[ii]
  #   w <- w[ii]
  #   wbar <- w/mean(w)
  #   S <- cumsum(wbar)
  #   pk <- (S + 1/2 - w/2)/(n + 1)
  #   approx(pk, x, p, method = "constant", f = 0, type = 2)$y
  # }
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
  posnext <- if (pos == n) pos else pos+1
  
  list(qlow = x[pos], qup = x[posnext],
       ilow = pos, iup = posnext,
       wlow = p - cumw[pos]/sum(w), # error?
       wup = cumw[posnext]/sum(w) - p)
}
