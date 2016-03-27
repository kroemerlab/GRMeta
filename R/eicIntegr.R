
integrOneEic<-function(tmpeic,lSamp=NULL,ivMint,eicParams,whichrt="rtcor",whichmz="mzcor"){
  
  bw=eicParams$bw
  bws=eicParams$nsmo1*bw
  bws2=eicParams$nsmo2*bw
  
  
  m=.GRconvEIC(tmpeic,whichrt=whichrt,bw=bws,delrt=bw/2)
  lrt=m[[2]];m=m[[1]]
  mz=.GRconvEIC(tmpeic,whichrt=whichrt,bw=bws2,delrt=bw/2,xnew=lrt)[[1]]
  
  if(is.null(lSamp)) lSamp=colnames(m)
  lsa=colnames(m)[colnames(m)%in%lSamp]
  if(length(lsa)==0) return(NULL)
  cm=m[,lsa,drop=F]
  cmz=mz[,lsa,drop=F]
  if(ncol(cm)>1) tabres=cbind(rt=lrt,.GRcompSeg(cm,cmz,eicParams,ivMint=ivMint))
  if(ncol(cm)==1) tabres=cbind(rt=lrt,.GRcompSegOne(cm,cmz,eicParams,ivMint=ivMint))
  
  ###########################
  apks=lapply(colnames(m),function(x) 
    .GRmsPeakSimple3(m[,x],span=eicParams$nspan,noise=mz[,x],snr.thresh=2))
  l=which(sapply(apks,nrow)>0)
  apks=do.call("rbind",lapply(l,function(x) data.frame(samp=colnames(m)[x],apks[[x]],stringsAsFactors=F)))
  rownames(apks)=NULL
  
  inloc=inlocr=inlocl=rep(0,nrow(apks))
  lseg=sort(unique(tabres[tabres[,4]!=0,2]))
  lseg=lseg[lseg>0]
  for(ij in lseg) inloc[apks[,2]%in%which(tabres[,2]==ij)]=ij
  for(ij in lseg) inlocl[apks[,3]%in%which(tabres[,2]==ij)]=ij
  for(ij in lseg) inlocr[apks[,4]%in%which(tabres[,2]==ij)]=ij
  apks=cbind(apks,pkreg=inloc,lefreg=inlocl,rireg=inlocr)
  
  lpks=which(tabres[,4]!=0)
  lseg=tabres[lpks,"Segment"]
  for(k in unique(lseg)) names(lpks)[lseg==k]=names(lseg)[lseg==k]=paste(k,1:sum(lseg==k),sep=".")
  finalpks=.GRgroupPks(apks,lpks,lseg)
  ###########################
  ## clean eic
  ieic=tmpeic$eic[1]
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
    pkstats[[i]]=data.frame(samp=isamp,Eic=ieic,Pk=finalpks$uuregs[i],rbind(c(RT.le=le,RT.ri=ri,Areas,HE,HEsm,MZs,Coda=Coda2)),stringsAsFactors=F)
    
  }
  neweic=do.call("rbind",neweic)
  pkstats=do.call("rbind",pkstats)
  
  list(Eic=neweic,SampStats=pkstats,Eicstats=NULL)
}
