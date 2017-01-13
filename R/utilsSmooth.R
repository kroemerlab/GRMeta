
.GRdoksmooth<-function(eic,bw){
  x=eic$x
  y=eic$y
  y2=(y!=0)*1
  if(!is.null(eic$cons)) lmiss=which(y==0 & eic$cons)
  else lmiss=which((y2[2:(length(y)-1)]-y2[1:(length(y)-2)]-y2[3:(length(y))])==-2)+1
  
  if(length(lmiss)>0){
    y2=c(y[1],y,y[length(y)])
    rmse=10^12
    y2[lmiss+1]=(y2[lmiss]+y2[lmiss+2])/2
    y2=y2[-c(1,length(y2))]
    del=10
    it=0
    while(del>10^-4){
      it=it+1
      kres=ksmooth(x, y2, kernel = "normal", bandwidth = bw,x.p=x)
      del=abs(rmse-mean((y2[lmiss]-kres$y[lmiss])^2))
      rmse=mean((y2[lmiss]-kres$y[lmiss])^2)
      y2[lmiss]=kres$y[lmiss]
    }
    y=y2
  }
  kres=ksmooth(x, y, kernel = "normal", bandwidth = bw)
  return(kres)
}

##### updated with over than zero
.GRdoksmooth2=function (eic, bw,minzero=0) 
{
  x = eic$x
  y = eic$y
  y2 = (y > minzero) * 1
  if (!is.null(eic$cons)) 
    lmiss = which(y <=minzero & eic$cons)
  else lmiss = which((y2[2:(length(y) - 1)] - y2[1:(length(y) - 
                                                      2)] - y2[3:(length(y))]) == -2) + 1
  if (length(lmiss) > 0) {
    y2 = c(y[1], y, y[length(y)])
    rmse = 10^12
    y2[lmiss + 1] = (y2[lmiss] + y2[lmiss + 2])/2
    y2 = y2[-c(1, length(y2))]
    del = 10
    it = 0
    while (del > 10^-4) {
      it = it + 1
      kres = ksmooth(x, y2, kernel = "normal", bandwidth = bw, 
                     x.p = x)
      del = abs(rmse - mean((y2[lmiss] - kres$y[lmiss])^2))
      rmse = mean((y2[lmiss] - kres$y[lmiss])^2)
      y2[lmiss] = kres$y[lmiss]
    }
    y = y2
  }
  kres = ksmooth(x, y, kernel = "normal", bandwidth = bw)
  return(kres)
}

########################################################################
##  baseline from ptw: 1.9-11
.GRasysm<-function (y, lambda = 1e+07, p = 0.001, eps = 1e-08, maxit = 25){
  z <- 0 * y
  w <- z + 1
  eps <- max(eps, eps * diff(range(y)))
  for (it in 1:maxit) {
    zold <- z
    z <- .GRwhit2(y, lambda, w)
    w <- p * (y > z) + (1 - p) * (y <= z)
    dz <- max(abs(z - zold))
    if (dz < eps) 
      break
  }
  if (dz >= eps) 
    warning("Function asysm did not reach convergence")
  return(z)
}

.GRwhit2<-function (y, lambda, w = rep(1, ny)){
  ny <- length(y)
  z <- d <- c <- e <- rep(0, length(y))
  .C("smooth2", w = as.double(w), y = as.double(y), z = as.double(z), 
     lamb = as.double(lambda), mm = as.integer(length(y)), 
     d = as.double(d), c = as.double(c), e = as.double(e), 
     PACKAGE = "GRMeta")$z
}
#################################
.GRbslrf<-function (x, y, span = 2/3, NoXP = NULL, maxit = c(2, 2), b = 3.5, 
          weight = NULL, Scale = function(r) median(abs(r))/0.6745, 
          delta = NULL, SORT = TRUE, DOT = FALSE, init = NULL) 
{
  if (!is.null(NoXP)) {
    if (NoXP < 3) 
      stop("NoXP is too small")
    NoXP <- as.integer(NoXP)
  }
  n <- length(x)
  if (setDelta <- is.null(delta)) {
    if (n <= 100) 
      delta <- 0
    else delta <- diff(range(x))/100
  }
  if (is.null(weight))  weight <- rep(1, n)
  
  if (is.null(init)) {
    rw <- rep(1, n)
    MAXIT <- c(maxit[1], sum(maxit)) + 1
  }
  else {
    if (length(init) != n) 
      stop("Must have same number of initial values as observations")
    Resid <- (y - init) * sqrt(weight)
    scale <- Scale(Resid)
    if (scale == 0) 
      status <- stop("could not compute scale of residuals for iter = 0")
    else {
      if (maxit[1] == 0) {
        u <- abs(Resid/(scale * b))
        rw <- ifelse(Resid < 0, 1, ifelse(u > 1, 0, 
                                          ((1 + u) * (1 - u))^2))
      }
      else {
        u <- abs(Resid/(scale * b))
        rw <- ifelse(u > 1, 0, ((1 + u) * (1 - u))^2)
      }
    }
    MAXIT <- c(maxit[1], sum(maxit))
  }
  if (SORT) {
    h <- order(x)
    x <- x[h]
    y <- y[h]
    weight <- weight[h]
    rw <- rw[h]
  }
  if (!is.loaded("lwreg")) {
    dyn.load("lwreg.dll")
    on.exit(dyn.unload("lwreg.dll"))
  }
  fit <- NULL
  for (iiter in 1:MAXIT[2]) {
    # if (DOT) {
    #   ok <- rw > 10^(-6)
    #   xx <- x[ok]
    #   wweight <- weight[ok] * rw[ok]
    #   yy <- y[ok]
    #   nn <- sum(ok)
    #   if (setDelta) {
    #     if (nn <= 100) 
    #       delta <- 0
    #     else delta <- diff(range(x))/100
    #   }
    #   qub <- if (is.null(NoXP)) 
    #     max(floor(span * nn), 2)
    #   else NoXP
    #   fit <- .C("lwreg", PACKAGE = "GRMeta", xx = as.double(xx), 
    #             as.double(yy), as.integer(nn), as.integer(qub), 
    #             as.double(delta), as.double(wweight), yfit = double(nn))
    # }
    # else {
      qub <- if (is.null(NoXP))
        max(floor(span * n), 2)
      else NoXP
      wweight <- weight * rw
      fit <- .C("lwreg", PACKAGE = "GRMeta", xx = as.double(x), 
                as.double(y), as.integer(n), as.integer(qub), 
                as.double(delta), as.double(wweight), yfit = double(n))
  #  }
    if (iiter < MAXIT[2]) {
      Resid <- (y - approx(x = fit$xx, y = fit$yfit, xout = x, 
                           rule = 2)$y) * sqrt(weight)
      scale <- Scale(Resid)
      if (scale == 0) 
        status <- stop(paste("could not compute scale of residuals for iter =", 
                             iiter))
      else {
        if (iiter < MAXIT[1]) {
          u <- abs(Resid/(scale * b))
          rw <- ifelse(Resid < 0, 1, ifelse(u > 1, 0, 
                                            ((1 + u) * (1 - u))^2))
        }
        if (iiter >= MAXIT[1]) {
          u <- abs(Resid/(scale * b))
          rw <- ifelse(u > 1, 0, ((1 + u) * (1 - u))^2)
        }
      }
    }
  }
  list(x = x, y = y, fit = approx(x = fit$xx, y = fit$yfit, 
                                  xout = x, rule = 2)$y, rw = rw, sigma = scale, span = span, 
       NoXP = NoXP, maxit = maxit, b = b, weight = weight, 
       Scale = Scale, delta = delta, SORT = SORT, DOT = DOT)
}
