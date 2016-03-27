.GRcompVals<-function(s,r,lenout=max(length(s),length(r)),maxsh=0){
  
  nlen=max(length(r),length(s))
  r=c(r,rep(0,lenout-length(r)))
  s=c(s,rep(0,lenout-length(s)))
  fftR <- fft(r)
  fftS <- fft(s)
  R <- fftR * Conj(fftS)
  R <- R/lenout
  rev = fft(R, inverse = TRUE)/sqrt(sum(r^2))/sqrt(sum(s^2))
  rxy=Re(rev)
  N=ceiling(lenout/2)
  rxy=rxy[c((N+1):lenout,1:(N))]
  if(maxsh>0){
    rxy=rxy[(N+1-min(maxsh,N)):(N+1+min(maxsh,N))]
    if(maxsh>N) rxy=c(rep(0,maxsh-N),rxy,rep(0,maxsh-N))
    rxy[is.na(rxy)]=0
  }
  rxy
}

.GRcompTICaddRT<-function(rtmat,intmat,rtrange=c(0,100),lref=colnames(rtmat)[grep("QC",colnames(rtmat))],drt=NULL,maxshift=NULL){ 
  
  if(is.null(drt)) drt=round(median(apply(mattic$rt,2,function(x) median(diff(x),na.rm=T))),4)
  lref=lref[lref%in%colnames(rtmat)]
  if(is.null(lref) | length(lref)<1) lref=sample(colnames(rtmat))[1:min(10,ncol(rtmat))]
  
  if(is.null(maxshift)) maxshift=20*drt
  rtgrid=round(rtmat/drt)
  lsc=max(floor(rtrange[1]/drt)-3,min(rtgrid,na.rm=T)):min(ceiling(rtrange[2]/drt)+3,max(rtgrid,na.rm=T))
  
  ## common grid
  intgrid=list()
  for(i in colnames(rtgrid)){
    iinit=tapply(intmat[,i],rtgrid[,i],sum,na.rm=T)
    iinit=iinit[match(lsc,as.numeric(names(iinit)))]
    l=which(is.na(iinit))
    if(length(l)>0)  iinit[l]=approx(lsc[-l],iinit[-l],lsc[l],yleft=iinit[-l][1],yright=rev(iinit[-l])[1])$y
    intgrid[[i]]=iinit
  }
  intgrid=do.call("cbind",intgrid)
  
  ################
  cat("\nComp shiftfrom: ",lref,"\n",sep=" ")
  lenout=2^ceiling(log2(nrow(intgrid)*2))
  maxsh=ceiling(maxshift/drt)
  if(maxsh%%2==0) maxsh=maxsh+1
  lsh=c(-maxsh:maxsh)*drt
  
  winM=corM=list()
  for(iref in lref){
    rec=apply(sqrt(intgrid),2,.GRcompVals,sqrt(intgrid)[,iref],lenout=lenout,maxsh=maxsh)
    rownames(rec)=lsh
    rec=rec[order(abs(lsh)),,drop=F]
    lx=as.numeric(rownames(rec))
    winM[[iref]]=apply(rec,2,function(x) lx[which.max(x)])
    corM[[iref]]=apply(rec,2,max,na.rm=T)
  }
  addrt=do.call("cbind",winM)
  corrrt=do.call("cbind",corM)
  return(list(AddRT=addrt,CorrRT=corrrt,IntGrid=intgrid))
}


.GRfiltreScan<-function(csc,idx=rep(1,length(csc)),alpha=5:10,perc=0.7){
  
  cmpmiss<-function(i) sum(i)/length(i)
  
  alpha=alpha[alpha<=min(max(alpha),length(csc))]
  isFeat=rep(FALSE,length(csc))
  for(add in alpha){
    z <- embed(csc,add)
    lrep=rowSums(z)/add
    l=which(lrep >=perc)
    l=l[which(idx[l]==idx[l+add-1])]
    for(i in 0:(add-1)) isFeat[l+i]=TRUE
  }
  return(isFeat)
}

.GRcodadw2<-function(x,lag=1) sum(diff(diff(x),lag=lag)^2)/sum(diff(x)^2)
.GRcodadw<-function(x) sum(diff(x)^2)/sum((x)^2)


