.GRchksign<-function(x,lim=10^-4){
  rex=re=rep(0,length(x))
  re[which(x>lim)]=1
  re[which(x< -lim)]=-1
  rex[which(diff(re)!=0)]=rex[which(diff(re)!=0)+1]=1
  rex[which(rex!=0)]=re[which(rex!=0)]+2
  rex
}

.GRatan=function(y) ifelse(y<0,1+atan(-y),atan(y))*180/pi

.GRgetArea<-function(x,y){
  n <- length(x)
  if(n==0) return(c(Sec=NA,Scan=NA))
  index <- order(x)
  dx <- diff(sort(x))
  z <- y[index]
  ys <- (z[1:(n - 1)] + z[2:n])/2
  return(c(Sec=sum(ys * dx),Scan=sum(ys)))
}

.GRmsExtrema<-function (x, span = 3){
  index1 <- .GRipeaks(x, span = span, strict = FALSE)
  index2 <- .GRipeaks(-x, span = span, strict = FALSE)
  index.max <- index1 & !index2
  index.min <- index2 & !index1
  iz=max(min(which(x>0))-1,1)
  if(!any(index.max[1:iz])) index.min[1:iz]=FALSE;
  index.min[iz]=TRUE;
  iz=min(length(x),max(which(x>0)+1))
  if(!any(index.max[iz:length(x)])) index.min[iz:length(x)]=FALSE
  index.min[iz]=TRUE;
  
  list(index.max = index.max, index.min = index.min)
}

.GRipeaks<-function (x, span = 3, strict = TRUE) 
{
  z <- embed(rev(as.vector(x)), dimension = span)
  z <- z[rev(seq(nrow(z))), ]
  s <- span%/%2
  v <- apply(z,1,which.max) == 1 + s
  z <- c(rep(FALSE, s), v)
  c(z[1:(length(z) - s)], rep(FALSE, span - 1))
}
