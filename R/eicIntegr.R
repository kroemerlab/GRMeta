.GRinfctInt<-function(tmp,m,isamp,lrt,le,ri,whichrt='whichrt',whichmz='whichmz'){
  lesm=which(lrt>=le & lrt<=ri)
  linpk=which(tmp$InPk & !is.na(tmp[,whichmz]))
  iapex=linpk[which.max(tmp$y[linpk])]
  MZs=c(tmp[iapex,whichmz],weighted.median(tmp[linpk,whichmz],tmp[linpk,"y"]))
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
  NCons=ifelse(is.null(tmp$NCons),NA,sum(tmp$NCons & tmp$InPk))
  reint=rbind(c(RT.le=le,RT.ri=ri,Areas,HE,HEsm,MZs,Coda=Coda2,NCons=NCons))
  return(reint)
}

integrOneEic<-function(tmpeic,lSamp=NULL,ivMint,eicParams,whichrt="rtcor",whichmz="mzcor",verbose=TRUE){
  
  
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
  
  if(all(tabres[,4]==0 | tabres[,2]==0)){
    if(verbose) cat("No master peaks in ",ieic,"\n")
    return(NULL)
  }
  ###########################
  apks=lapply(colnames(m),function(x) 
    .GRmsPeakSimple3(m[,x],span=eicParams$nspan,noise=mz[,x],snr.thresh=2))
  l=which(sapply(apks,nrow)>0)
  apks=do.call("rbind",lapply(l,function(x) data.frame(samp=colnames(m)[x],apks[[x]],stringsAsFactors=F)))
  rownames(apks)=NULL
  if(nrow(apks)==0){
    if(verbose) cat("No peaks in sample",ieic,"\n")
    return(NULL)
  }
  
  finalpks=.GRgroupPks(apks,tabres)
  finalpks$IsNA=is.na(finalpks$Pk)
  if(any(is.na(finalpks$Pk))){
    l=which(is.na(finalpks$Pk))
    finalpks$Pk[l]=max(finalpks$Pk,na.rm=T)+(1:length(l))
  }
  if(verbose){
    if(is.null(finalpks$Pk)) cat("No peaks in ",ieic,"\n")
    if(any(finalpks$IsNA)) cat("NAs in ",ieic,"\n")
  }
  finalpks=finalpks[which(finalpks$Pk!=0),,drop=F]
  if(nrow(finalpks)==0) return(NULL)
  finalpks0=finalpks
  ###########################
  ## reduce final list
  scl=paste(finalpks0$samp,finalpks0$Pk,sep=";;")
  finalpks=data.frame(samp=tapply(finalpks0$samp,scl,unique),
                      tick.loc=tapply(1:nrow(finalpks0),scl,function(x) finalpks0$tick.loc[x][which.max(finalpks0$SNR[x])]),
                      tick.left=tapply(finalpks0$tick.left,scl,min),
                      tick.right=tapply(finalpks0$tick.right,scl,max),
                      SNR=tapply(finalpks0$SNR,scl,max),
                      IsNA=tapply(finalpks0$IsNA,scl,all),
  Pk=tapply(finalpks0$Pk,scl,unique),stringsAsFactors=FALSE)
  for(i in names(finalpks)) finalpks[,i]=as.vector(finalpks[,i])
  
  
  ###########################
  ## get final pks stats
  finalpks$MZ=finalpks$MZ2=finalpks$RT=NA
  for(i in 1:nrow(finalpks)){
    le=lrt[finalpks$tick.left[i]]
    ri=lrt[finalpks$tick.right[i]]
    isamp=finalpks$samp[i]
    tmp=tmpeic[tmpeic[,whichrt]>=(le) & tmpeic[,whichrt]<=(ri) & tmpeic$samp==isamp,]
    iapex=which.max(tmp$y)
    finalpks$RT[i]=tmp[iapex,whichrt]
    finalpks$MZ[i]=tmp[iapex,whichmz]
    finalpks$MZ2[i]=weighted.median(tmp[,whichmz],tmp[,"y"])
  }
  if(max(finalpks$Pk)>1){
    ############################
    ## reorder based on RT/MZ
    RTs=round(tapply(finalpks$RT,finalpks$Pk,median),4)
    MZs=tapply(finalpks$MZ,finalpks$Pk,median)
    lso=names(RTs)[order(RTs,MZs)]
    finalpks$Pk=(1:length(lso))[match(finalpks$Pk,lso)]
    
    lcl=as.list(1:length(lso))
    oRTs=sapply(lcl,function(x) median(finalpks$RT[finalpks$Pk%in%x]))
    oMZs=sapply(lcl,function(x) median(finalpks$MZ[finalpks$Pk%in%x]))
    itab=cbind(diff(oRTs),abs(.GRcompdppm(oMZs,F)))
    
    deltaRT=ifelse(is.null(eicParams$dRT),eicParams$nspan*eicParams$bw*3,eicParams$dRT)
    deltaPPM=ifelse(is.null(eicParams$dPPM),11,eicParams$dPPM)
    while(any((itab[,1]<deltaRT & itab[,2]<deltaPPM))){
      l=which(itab[,1]<deltaRT & itab[,2]<deltaPPM)
      l=l[which.min(itab[l,1])]
      lcl[[l+1]]=unlist(lcl[l:(l+1)])
      lcl=lcl[-l]
      if(length(lcl)==1) break
      oRTs=sapply(lcl,function(x) median(finalpks$RT[finalpks$Pk%in%x]))
      oMZs=sapply(lcl,function(x) median(finalpks$MZ[finalpks$Pk%in%x]))
      itab=cbind(diff(oRTs),abs(.GRcompdppm(oMZs,F)))
    }
    oldPk=finalpks$Pk
    for(i in 1:length(lcl)) finalpks$Pk[oldPk%in%lcl[[i]]]=i
    if(any(finalpks$IsNA)){
      lpks2rm=unique(finalpks$Pk[finalpks$IsNA])
      lpks2rm=lpks2rm[!lpks2rm%in%unique(finalpks$Pk[!finalpks$IsNA])]
      if(length(lpks2rm)>0) finalpks=finalpks[!finalpks$Pk%in%lpks2rm,,drop=F]
    }
    RTs=round(tapply(finalpks$RT,finalpks$Pk,median),4)
    MZs=tapply(finalpks$MZ,finalpks$Pk,median)
    lso=names(RTs)[order(RTs,MZs)]
    finalpks$Pk=(1:length(lso))[match(finalpks$Pk,lso)]
    finalpks1=finalpks
    scl=paste(finalpks1$samp,finalpks1$Pk,sep=";;")
    finalpks=data.frame(samp=tapply(finalpks1$samp,scl,unique),
                        tick.loc=tapply(1:nrow(finalpks1),scl,function(x) finalpks1$tick.loc[x][which.max(finalpks1$SNR[x])]),
                        tick.left=tapply(finalpks1$tick.left,scl,min),
                        tick.right=tapply(finalpks1$tick.right,scl,max),
                        SNR=tapply(finalpks1$SNR,scl,max),
                        IsNA=tapply(finalpks1$IsNA,scl,all),
                        Pk=tapply(finalpks1$Pk,scl,unique),stringsAsFactors=FALSE)
    for(i in names(finalpks)) finalpks[,i]=as.vector(finalpks[,i])
    
  }
  
  ###########################
  ## clean eic
  bws2=eicParams$nsmo2*eicParams$bw
  neweic=pkstats=list()
  for(i in 1:nrow(finalpks)){
    le=lrt[finalpks$tick.left[i]]
    ri=lrt[finalpks$tick.right[i]]
    isamp=finalpks$samp[i]
    tmp=tmpeic[tmpeic[,whichrt]>=(le-bws2/2) & tmpeic[,whichrt]<=(ri+bws2/2) & tmpeic$samp==isamp,]
    tmp$NCons=.GRfiltreScan(tmp$y>=ivMint,alpha=eicParams$LSc,perc=eicParams$Perc)
    tmp$InPk=(tmp[,whichrt]>=le & tmp[,whichrt]<=ri)
    tmp$Pk=finalpks$Pk[i]
    neweic[[i]]=tmp
    ### lrt,le,ri,tmp,m
    resint=.GRinfctInt(tmp,m,isamp,lrt,le,ri,whichrt,whichmz)
    pkstats[[i]]=data.frame(samp=isamp,eic=ieic,Pk=finalpks$Pk[i],resint,stringsAsFactors=F)
    
  }
  neweic=do.call("rbind",neweic)
  pkstats=do.call("rbind",pkstats)
  if(length(unique(pkstats$Pk))==1) return(list(Eic=neweic,SampStats=pkstats,Eicstats=NULL))
  # neweic=allre[[ieic]]$Eic ; pkstats=allre[[ieic]]$SampStats
  
  list(Eic=neweic,SampStats=pkstats,Eicstats=NULL)
}

####################
.GRintegrOneEicGrpCl<-function(tabeic,lSamp=NULL,eicParams,whichrt="rtcor",whichmz="mzcor",ncl=4){
 
  .inGRintegrOneEicGrpCl<-function(igrpeic,tabeic,lSamp,eicParams,whichrt,whichmz){
     itabeic=tabeic[tabeic$GrpEic==igrpeic,]
     tt=integrOneEicGrp(itabeic,lSamp,eicParams,whichrt,whichmz,save=TRUE,verbose=FALSE)
  }
  
   lgrpeic=unique(tabeic$GrpEic)
   lfiles=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
               ieicgrp,
               ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
   
   lgrpeic=lgrpeic[file.exists(lfiles)]
   tabeic=tabeic[tabeic$GrpEic%in%lgrpeic,]
   d0=proc.time()[3]
   cat("Started at ",date(),sep="")
     require("snowfall")
     ncl=max(1,min(ncl,parallel:::detectCores()))
     cat(" on ",ncl," processors\n",sep="")
     sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
     sfExport( ".inGRintegrOneEicGrpCl", local=TRUE )
     sfLibrary(GRMeta)
     allr=sfClusterApplyLB(lgrpeic,.inGRintegrOneEicGrpCl,tabeic,lSamp,eicParams,whichrt,whichmz)
     sfStop()
   d1=proc.time()[3]
   cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs to process ",length(lgrpeic)," EIC groups\n",sep="")
   
  
}
  

integrOneEicGrp<-function(tabeic,lSamp=NULL,eicParams,whichrt="rtcor",whichmz="mzcor",save=TRUE,verbose=TRUE){
  
  ieicgrp=tabeic$GrpEic[1]
  tabeic=tabeic[tabeic$GrpEic==ieicgrp,]
  
  ifile=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
              ieicgrp,
              ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
  
  load(ifile)
  
  dfeic=dfeic[dfeic$eic%in%tabeic$Id,]
  rownames(dfeic)=NULL
  eicdef=attr(dfeic,'eic')
  attr(dfeic,'eic')=eicdef=eicdef[eicdef$Id%in%dfeic$eic,]
  leics=unique(dfeic$eic)
  if(verbose) cat("Found ",length(leics)," in ",ieicgrp,"\n",sep="")
  allre=list()
  for(ieic in leics){
    itmpeic=dfeic[dfeic$eic==ieic,]
    ivMint=tabeic$Bl[tabeic$Id==ieic]
    ivMint=ifelse(is.null(ivMint),eicParams$Mint,ivMint)
    if(nrow(itmpeic)>0){
      if(verbose) cat(" * ",ieic,": ",sep="")
      allre[[ieic]]=integrOneEic(itmpeic,lSamp=lSamp,ivMint=ivMint,eicParams=eicParams,
                            whichrt=whichrt,whichmz=whichmz,verbose=verbose)
      if(verbose) 
      cat(ifelse(is.null(allre[[ieic]]),' 0 ',length(unique(allre[[ieic]]$SampStats$Pk)))," peaks\n",sep="")
    }
  }
  
  allre=allre[which(!sapply(allre,is.null))]
  cat("\n")
  if(length(allre)==0) return(NULL)
  
  dfeic=do.call("rbind",lapply(allre,function(x) x$Eic))
  dfeic$eic=paste(dfeic$eic,dfeic$Pk,sep="-")
  dfeic=dfeic[,names(dfeic)!="Pk"]
  
  sstats=do.call("rbind",lapply(allre,function(x) x$SampStats))
  sstats$eic0=sstats$eic
  sstats$eic=paste(sstats$eic,sstats$Pk,sep="-")

  rm(list=c("allre"))
  ############## Combine
  lsids=sort(unique(sstats$samp))
  lsids1=lsids[lsids%in%lSamp]
  lsids2=lsids[!lsids%in%lSamp]
  
  leics=sort(unique(sstats$eic))
  sstats=sstats[order(factor(sstats$eic,levels=leics),factor(sstats$samp,lsids)),,drop=F]
  ldats=names(sstats)[!names(sstats)%in%c("samp","eic","eic0","Pk")]
  alldat=vector('list',length(ldats));names(alldat)=ldats
  allstats=alldat
  convm=do.call("cbind",tapply(1:nrow(sstats),sstats$eic,function(x) x[match(lsids,sstats$samp[x])]))
  rownames(convm)=lsids
  lnna=which(!is.na(convm))
  for(idat in ldats){
    m=matrix(NA,nrow=length(lsids),ncol=length(leics),dimnames=dimnames(convm))
    m[lnna]=sstats[convm[lnna],idat]
    ist=rep(NA,length(leics))
    if(length(lsids1)>0) ist=apply(m[lsids1,,drop=F],2,median,na.rm=T)[leics]
    if(length(lsids2)>0 & any(is.na(ist))) ist[which(is.na(ist))]=apply(m[lsids1,which(is.na(ist)),drop=F],2,median,na.rm=T)
    allstats[[idat]]=ist
    alldat[[idat]]=m
    
  }
  allstats=t(do.call("rbind",allstats))
  allstats=data.frame(eic=rownames(allstats),allstats,stringsAsFactors=FALSE)
  
  eicstats=data.frame(PkId=rownames(allstats),eicori=as.vector(tapply(sstats$eic0,sstats$eic,unique)[allstats$eic]),
   pkori=as.vector(tapply(sstats$Pk,sstats$eic,unique)[allstats$eic]),stringsAsFactors=FALSE)
  eicstats=cbind(eicstats[,-2],tabeic[match(eicstats$eicori,tabeic$Id),])
  rownames(eicstats)=eicstats$PkId
  
  rm('sstats')
  if(save){
    ofile=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
                ieicgrp,
                ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),"-clean.rda",sep="")
    save(file=ofile,dfeic,alldat,allstats,eicstats)
  }
  
  invisible(list(Eic=dfeic,Dat=alldat,Stats=allstats,EicInfos=eicstats))
  
}
