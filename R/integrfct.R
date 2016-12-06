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

.GRipeaks<-function (x, span = 3, strict = TRUE) 
{
  z <- embed(rev(as.vector(x)), dimension = span)
  z <- z[rev(seq(nrow(z))), ]
  s <- span%/%2
  v <- apply(z,1,which.max) == 1 + s
  z <- c(rep(FALSE, s), v)
  c(z[1:(length(z) - s)], rep(FALSE, span - 1))
}

.GRmsPeakSimple<-function (x, y, noise.local = NULL, span = 3, snr.thresh = 2) {
  
  index <- .GRmsExtrema(y, span = span)
  nvar <- length(x)
  
  index.min = index$index.min; index.max = index$index.max
  
  imax <- which(index.max)
  snr=y[imax]
  if(!is.null(noise.local)) snr=snr/noise.local[imax]
  good.snr <- (snr > snr.thresh)
  snr <- snr[good.snr]
  tick.loc <- imax[good.snr]
  if ((npeak = length(tick.loc)) == 0) return(data.frame())
  tick.left <- tick.right <- rep(1, length(tick.loc))
  for (i in 1:length(tick.loc)) {
    #for (j in tick.loc[i]:1) 
    tick.left[i] <-max(which(index.min[1:tick.loc[i]]))
    #for (j in tick.loc[i]:length(x)) 
    tick.right[i] <- min(which(index.min[tick.loc[i]:length(x)]))+tick.loc[i]-1
  }
  if (tick.right[length(tick.loc)] == 1) tick.right[length(tick.loc)] <- length(x)
  if (tick.left[1] == 1) tick.left[1] <- which.min(y[1:tick.loc[1]]!=0)
  
  for(ix in 1:length(tick.loc)){
    l=tick.left[ix]:tick.right[ix]
    tick.left[ix]=max(tick.left[ix],min(l[which(y[l]>0)])-1)
    tick.right[ix]=min(tick.right[ix],max(l[which(y[l]>0)])+1)
  }
  ###  check duplicated left -> use the highest peak
  lpks=lapply(unique(tick.left),function(x) which(tick.left==x))
  pks=do.call("rbind",lapply(lpks,function(ix){
    if(length(ix)>1) return(c(tick.loc=tick.loc[ix[which.max(y[tick.loc[ix]])]],tick.left=min(tick.left[ix]), tick.right=max(tick.right[ix])))
    c(tick.loc=tick.loc[ix],tick.left=tick.left[ix], tick.right=tick.right[ix])
  }))
  pks=data.frame(pks,tick.span = pks[,3]-pks[,2] + 1, snr = y[pks[,1]], intensity = y[pks[,1]],
                 mass.loc = x[pks[,1]], mass.left=x[pks[,2]], mass.right=x[pks[,3]], mass.span=x[pks[,3]]-x[pks[,2]])
  ### compute clusters
  icl=0;lcl=l=c()
  for(ix in 1:nrow(pks)){
    if(any((pks[ix,3]:pks[ix,2])%in%l)){
      l=c(l,pks[ix,3]:pks[ix,2])
    }
    else{
      l=pks[ix,3]:pks[ix,2]
      icl=icl+1      
    }
    lcl=c(lcl,icl)    
  }
  pks$cl=lcl
  
  return(pks)
  
}

.GRconvEIC<-function(tmp,whichrt="rt",delrt=0.005,bw=0.025,xnew=NULL){
  if(is.null(xnew))
    xnew=seq(floor(min(tmp[,whichrt]*10))/10,ceiling(max(tmp[,whichrt]*10))/10,delrt)
  amat=NULL
  amz=NULL
  lsamp=sort(unique(tmp[,"samp"]))
  for(k in lsamp){
    l=which(tmp[,"samp"]==k)
    xl=tmp[l,whichrt]
    yl=tmp[l,"y"]
    drt=quantile(diff(xl),c(.25,.50));drt=unique(drt[drt>0])[1]
    if(is.na(drt)) drt=quantile(diff(xnew),.25)
    padl=max(2,ceiling((min(xl)-min(xnew))/drt))
    padr=max(2,floor((max(xnew)-max(xl))/drt)+1)
    yl=c(rep(0,padl),yl,rep(0,padr))
    xl=c(min(xl)-(padl:1)*drt,xl,(1:padr)*drt+max(xl))
    kres=.GRdoksmooth(list(x=xl,y=yl,cons=NULL),bw)
    amat=cbind(amat,approx(kres$x,kres$y,xnew)$y)
  }
  amat[is.na(amat)]=0
  colnames(amat)=lsamp
  rownames(amat)=xnew
  return(list(y=amat,x=xnew))
}
