
integrOneEic<-function(tmpeic,lSamp=NULL,ivMint,eicParams,whichrt="rtcor",whichmz="mzcor"){
  
  
  ieic=tmpeic$eic[1]
  m=.GRconvEIC(tmpeic,whichrt=whichrt,bw=eicParams$nsmo1*eicParams$bw,delrt=eicParams$bw/2)
  lrt=m[[2]];m=m[[1]]
  mz=.GRconvEIC(tmpeic,whichrt=whichrt,bw=eicParams$nsmo2*eicParams$bw,delrt=eicParams$bw/2,xnew=lrt)[[1]]
  
  if(is.null(lSamp)) lSamp=colnames(m)
  lsa=colnames(m)[colnames(m)%in%lSamp]
  if(length(lsa)==0) return(NULL)
  cm=m[,lsa,drop=F]
  cmz=mz[,lsa,drop=F]
  if(ncol(cm)>1) tabres=cbind(rt=lrt,.GRcompSeg(cm,cmz,eicParams,ivMint=ivMint))
  if(ncol(cm)==1) tabres=cbind(rt=lrt,.GRcompSegOne(cm,cmz,eicParams,ivMint=ivMint))
  
  ###########################
  apks=lapply(colnames(m),function(x) 
    .GRmsPeakSimple3(m[,x],span=eicParams$nspan,noise=mz[,x],snr.thresh=1.1))
  l=which(sapply(apks,nrow)>0)
  apks=do.call("rbind",lapply(l,function(x) data.frame(samp=colnames(m)[x],apks[[x]],stringsAsFactors=F)))
  rownames(apks)=NULL
  
  finalpks=.GRgroupPks(apks,tabres)
  if(is.null(finalpks$Pk)) cat("No peaks in ",ieic,"\n")
  if(any(is.na(finalpks$Pk))) cat("NAs in ",ieic,"\n")
  finalpks=finalpks[which(finalpks$Pk>0),,drop=F]

  ###########################
  ## clean eic
  bws2=eicParams$nsmo2*eicParams$bw
  neweic=pkstats=list()
  for(i in 1:nrow(finalpks)){
    le=lrt[finalpks$tick.left[i]]
    ri=lrt[finalpks$tick.right[i]]
    isamp=finalpks$samp[i]
    tmp=tmpeic[tmpeic[,whichrt]>=(le-bws2/2) & tmpeic[,whichrt]<=(ri+bws2/2) & tmpeic$samp==isamp,]
    tmp$InPk=(tmp[,whichrt]>=le & tmp[,whichrt]<=ri)
    tmp$Pk=finalpks$uuregs[i]
    neweic[[i]]=tmp
    ###
    lesm=which(lrt>=le & lrt<=ri)
    linpk=which(tmp$InPk & !is.na(tmp[,whichmz]))
    iapex=linpk[which.max(tmp$y[linpk])]
    MZs=c(tmp[iapex,whichmz],weighted.median(tmp[linpk,whichmz],tmp[linpk,whichmz]))
    names(MZs)=c("MZap","MZwm")
    Areas=rep(.GRgetArea(lrt[lesm],m[lesm,isamp]),2)
    if(sum(tmp$InPk)>2) Areas[1:2]=.GRgetArea(tmp[tmp$InPk,whichrt],tmp[tmp$InPk,"y"])
    names(Areas)=c("Area.rt","Area.scan","Area.rtsm","Area.scansm")
    iapexsm=lesm[which.max(m[lesm,isamp])]
    HEsm=c(RT.apsm=lrt[iapexsm],HE.apsm=m[iapexsm,isamp])
    linpk=which(tmp$InPk)
    iapex=linpk[which.max(tmp$y[linpk])]
    HE=c(RT.ap=tmp[iapex,whichrt],HE.ap=tmp[iapex,"y"])
    Coda2=.GRcodadw2(tmp[tmp$InPk,'y'])
    Coda2[is.na(Coda2)]=3
    pkstats[[i]]=data.frame(samp=isamp,eic=ieic,Pk=finalpks$Pk[i],rbind(c(RT.le=le,RT.ri=ri,Areas,HE,HEsm,MZs,Coda=Coda2)),stringsAsFactors=F)
    
  }
  neweic=do.call("rbind",neweic)
  pkstats=do.call("rbind",pkstats)
  
  list(Eic=neweic,SampStats=pkstats,Eicstats=NULL)
}
