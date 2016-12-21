
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