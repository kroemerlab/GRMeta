
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

### collapse the eic attr
.GRgetEicInfos<-function(eicParams,eicloc=NULL){
  
  if(!is.null(eicloc)) if(is.data.frame(eicloc)){
  fileeics=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
                 unique(eicloc$GrpEic),
                 ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
  names(fileeics)=unique(eicloc$GrpEic)
  }
  
  if(is.null(eicloc))
    fileeics=list.files(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),pattern=paste(eicParams$addEic,".rda",sep=""),full.names=TRUE)
  
  fileeics=fileeics[file.exists(fileeics)]
  if(length(fileeics)==0) stop('No file found!')
  
  allr=list()
  for(k in 1:length(fileeics)){
    cat(ifelse(k%%10==0,"X","."))
    load(fileeics[k])
    allr[[k]]=attr(dfeic,"eic")
  }
  newtabeic=do.call("rbind",allr)
#  if(is.null(eicloc)){
    grp=factor(newtabeic$GrpEic,names(sort(tapply(newtabeic$mzap,newtabeic$GrpEic,mean,na.rm=T))))
    newtabeic=newtabeic[order(grp,newtabeic$mzap),]
 # }
  newtabeic
  
}

.GRestimBl<-function(tabeic,blSid,eicParams,ncl=1){
  
  .infctBlfct<-function(ifile,blSid=NULL,minScan=11,quantBl=0.5){
    load(ifile)
    dfeic=dfeic[dfeic$samp%in%blSid,]
    leics=attr(dfeic,"eic")[,"Id"]
    matbl=matrix(NA,nrow=length(leics),ncol=length(blSid),dimnames=list(leics,blSid))
    if("ineic" %in% names(dfeic)) dfeic=dfeic[which(dfeic$ineic==1),]
    if(nrow(dfeic)==0) return(matbl)
    eicbl=paste(dfeic$samp,dfeic$eic)
    l2use=names(which(table(eicbl)>=minScan))
    if(length(l2use)==0) return(matbl)
    dfeic=dfeic[which(eicbl%in%l2use & dfeic$y>0 & !is.na(dfeic$y)),]
    if(nrow(dfeic)==0) return(matbl)
    blsamp=tapply(1:nrow(dfeic),dfeic$samp,function(x){
      x=tapply(dfeic$y[x],dfeic$eic[x],quantile,quantBl)
      x=x[match(rownames(matbl),names(x))]
      names(x)=rownames(matbl)
      x
      })
    for(i in names(blsamp))    matbl[names(blsamp[[i]]),i]=blsamp[[i]]
    
    
 #   print("OKKK")
    rm(list=c("dfeic"))
    return(matbl)
  }
  
  
  fileeics=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
                 unique(tabeic$GrpEic),
                 ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
  names(fileeics)=unique(tabeic$GrpEic)
  fileeics=fileeics[file.exists(fileeics)]
  print(fileeics)
  if(length(fileeics)==0) stop('No file found!')
  
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  
  if(ncl==1 | length(fileeics)==1){
    cat(" on 1 processor\n",sep="")
    allr=list()
    for(k in 1:length(fileeics)){
      cx=fileeics[k]
      cat(cx," ")
      allr[[k]]=.infctBlfct(cx,blSid,eicParams$minScanBl,eicParams$quantBl)
    }
  }
  if(ncl>1 & length(fileeics)>1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfExport( ".infctBlfct", local=TRUE )
    allr=sfClusterApplyLB(fileeics,.infctBlfct,blSid,eicParams$minScanBl,eicParams$quantBl)
    sfStop()
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs \n",sep="")
  
  allr=allr[which(!sapply(allr,is.null))]
  cat("\n")
  if(length(allr)==0) return(NULL)
  blres=do.call("rbind",allr)
  blres=blres[rownames(blres)%in%tabeic$Id,,drop=F]
  blres=blres[match(tabeic$Id,rownames(blres)),,drop=F]
  rownames(blres)=tabeic$Id
  return(blres)
}


.GRcompEICst<-function(tabeic,eicParams,lSids,ncl=1){
  
  .incompEICst<-function(cfi,blvect=NULL,eicParams){
    
    load(cfi)  
    leics=attr(dfeic,"eic")$Id
    if(is.null(blvect)) blvect=eicParams$Mint
    if(length(blvect==1)){
      iblvect=rep(blvect,length(leics))
      names(iblvect)=leics
    } else iblvect=blvect[names(blvect)%in%leics]
    if(any(is.na(iblvect))) iblvect[is.na(iblvect)]=eicParams$Mint
    
    dfeic$csc=as.vector(dfeic$y>=iblvect[dfeic$eic])
    Ncons=Coda2=HMax=list()
    lusamp=unique(dfeic$samp)
    if("ineic" %in% names(dfeic)) unique(dfeic$samp[which(dfeic$ineic==1)])
    for(isamp in lusamp){
      idfeic=dfeic[which(dfeic$samp==isamp),]
      if("ineic" %in% names(idfeic)) idfeic=idfeic[which(idfeic$ineic==1),]
      re=.GRfiltreScan(idfeic$csc,idfeic$eic,alpha=eicParams$LSc,perc=eicParams$Perc) 
      Ncons[[isamp]]=tapply(re,idfeic$eic,sum)[leics]
      Coda2[[isamp]]=tapply(idfeic$y,idfeic$eic,.GRcodadw2)[leics]
      HMax[[isamp]]=tapply(idfeic$y,idfeic$eic,function(x) floor(max(x,na.rm=T)))[leics]
    }
    Ncons=do.call("rbind",Ncons)
    Coda2=do.call("rbind",Coda2)
    HMax=do.call("rbind",HMax)
    colnames(HMax)=colnames(Coda2)=colnames(Ncons)=leics
    return(list(Ncons=Ncons,Coda2=Coda2,HMax=HMax))
  }
  
  fileeics=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
                 unique(tabeic$GrpEic),
                 ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
  names(fileeics)=unique(tabeic$GrpEic)
  fileeics=fileeics[file.exists(fileeics)]
  if(length(fileeics)==0) stop('No file found!')
  
  if("Bl"%in%names(tabeic)) blvect=tabeic$Bl else blvect=rep(eicParams$Mint,nrow(tabeic))
  names(blvect)=tabeic$Id
  
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  
  if(ncl==1 | length(fileeics)==1){
    cat(" on 1 processor\n",sep="")
    allr=list()
    for(k in 1:length(fileeics)){
      cx=fileeics[k]
      cat(ifelse(k%%10==0,"X","."))
      allr[[k]]=.incompEICst(cx,blvect,eicParams)
    }
  }
  if(ncl>1 & length(fileeics)>1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfExport(".incompEICst", local=TRUE )
             sfLibrary(GRMeta)
             
#     sfExport( ".GRfiltreScan", local=TRUE )
#     sfExport( ".GRcodadw2", local=TRUE )
    allr=sfClusterApplyLB(fileeics,.incompEICst,blvect,eicParams)
    sfStop()
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs \n",sep="")
  
  allr=allr[which(!sapply(allr,is.null))]
  cat("\n")
  if(length(allr)==0) return(NULL)
  
  if(missing(lSids)) lSids=sort(unique(unlist(lapply(allr,function(x) rownames(x[[1]])))))
  allr=lapply(allr,function(x) lapply(x,function(y) y[match(lSids,rownames(y)),,drop=F]))
  
  lmats=unique(unlist(lapply(allr,names)))
  allmat=list()
 for(imat in lmats){
  tmpm=lapply(allr,function(x) x[[imat]])
  tmpm=tmpm[which(!sapply(tmpm,is.null))]
  tmpm=do.call("cbind",tmpm)
  tmpm=tmpm[,match(tabeic$Id,colnames(tmpm)),drop=F]
  dimnames(tmpm)=list(lSids,tabeic$Id)
  allmat[[imat]]=tmpm
 }
  rm(list="allr")
  return(allmat)
  
}




