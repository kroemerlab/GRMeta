
parseOneSampEIC<-function(mzfi,tabeic,outfile=NULL,npad=3,stepmz=1/1000,mzdata=FALSE,verbose=TRUE,serial=TRUE,chunk=200){
  
  
  d0=proc.time()[3]
  if(!file.exists(mzfi)) return(NULL)
  if(!"shrt"%in%names(tabeic)) tabeic$shrt=0
  if(!"shppm"%in%names(tabeic)) tabeic$shppm=0
  ##########################################
  ## load data
  if(mzdata){
    aa=openMSfile(mzfi,backend="Ramp",verb=T)
    rts=header(aa)[,"retentionTime"]/60
    names(rts)=header(aa)[,1]
    pl <- mzR::peaks(aa)
    mzR::close(aa)
    m=do.call("rbind",lapply(1:length(pl),function(x) cbind(Scan=rep(x,nrow(pl[[x]])),pl[[x]])))
    rm("pl")
    m=m[order(m[,2]),]
    m=cbind(m,rts[m[,1]])
    colnames(m)=c("scan","mz","y","rt")
  }else load(mzfi)
  
  whichmz='mz';whichrt='rt'
  
  tabeic=tabeic[order(tabeic[,'mzmin']),]
  
  
  padrt=ifelse(npad<0,0,npad*median(diff(rts))*(npad+1))
  
  ##########################################
  ## filter the data file: wide mz
  lcode=round(m[,whichmz]/stepmz)
  m=m[order(lcode,m[,whichrt]),]
  lcode=round(m[,whichmz]/stepmz)
  if(verbose) cat("Reduction from ",nrow(m)," ",sep="")
  
  llcode=split(1:nrow(tabeic), ceiling(seq_along(1:nrow(tabeic))/chunk))
  amred=list()
  for(ilins in 1:length(llcode)){
    if(verbose) cat(ifelse(ilins%%10==0,'*','.'))
    lins=llcode[[ilins]]
    itab=tabeic[lins,]
    rtmzrange=cbind(floor((itab[,'mzmin']*(1-itab[,'shppm']*10^-6)-2*stepmz)/stepmz),
                    ceiling((itab[,'mzmax']*(1+itab[,'shppm']*10^-6)+2*stepmz)/stepmz))
    l=which(lcode>=min(rtmzrange) & lcode<=max(rtmzrange))
    v=lcode[l]
    mred=m[l,]
#    system.time(rangmz<-apply(rtmzrange,1,function(x) xcms:::findRange(v,x[1],x[2])))
    #ucode=apply(t(rangmz),1,function(x) x[1]:x[2])
    system.time(ucode<-apply(rtmzrange,1,function(x) which(v>=x[1] & v <=x[2])))
     mred=mred[unlist(ucode),]
    mred=cbind(mred,eic=unlist(lapply(1:length(ucode),function(x) rep(x,length(ucode[[x]])))))
    mred=mred[which((mred[,whichrt]+itab[mred[,5],"shrt"]-itab[mred[,5],"rtmin"]+padrt)>=0),,drop=FALSE]
    mred=mred[which((mred[,whichrt]+itab[mred[,5],"shrt"]-itab[mred[,5],"rtmax"]-padrt)<=0),,drop=FALSE]
    mred=mred[which((mred[,whichmz]*(1+itab[mred[,5],"shppm"]*10^-6)-itab[mred[,5],"mzmin"])>=0),,drop=FALSE]
    mred=mred[which((mred[,whichmz]*(1+itab[mred[,5],"shppm"]*10^-6)-itab[mred[,5],"mzmax"])<=0),,drop=FALSE]
    mred[,"eic"]=lins[mred[,"eic"]]
    amred[[ilins]]=mred
    rm(list="mred")
  }
  rm(list=c("m",'lcode','llcode'))
  mred=do.call("rbind",amred)
  rm(list=c('amred'))
  mred<-mred[order(mred[,"eic"],mred[,"scan"]),]
  mred[which(is.na(mred[,'y'])),'y']=0
  mred[,'rt']=rts[as.character(mred[,'scan'])]
  if(verbose) cat(" plus ",nrow(mred),"\n",sep="")
  
  rownames(mred)=NULL  
  attr(mred,"eic")<-tabeic
  attr(mred,"file")<-mzfi
  attr(mred,"rts")<-rts
  
  d1=proc.time()[3]
  n=table(mred[,'eic'])
  top=paste(mzfi," :"," -> ",round(d1-d0,1)," secs :",length(n),"/",nrow(tabeic),": ",round(mean(n),3)," scan "," (",round(100*length(n)/nrow(tabeic),1),"%) EICs founds.",sep="")
  if(!is.null(outfile)) save(file=outfile,mred)
  if(!serial) return(invisible(list(top,d1-d0)))
  if(verbose)   cat(top,sep="\n")
  return(invisible(mred))
  
}


.GRparseSampEICCl<-function(ix,matfile,corrt,tabeic,npad,stepmz,mzdata,chunk){
  tmpeic=tabeic
  if("shppm"%in%names(matfile)) tmpeic$shppm=matfile[ix,]$shppm else tmpeic$shppm=0
  if("shrt"%in%names(matfile)) tmpeic$shrt=matfile[ix,]$shrt else tmpeic$shrt=0
  if(!is.null(corrt)) if(ix%in%rownames(corrt))
    tmpeic$shrt=tmpeic$shrt+approx(as.numeric(colnames(corrt)),corrt[ix,],tmpeic$rtap)$y
  
 parseOneSampEIC(matfile[ix,]$In,outfile=matfile[ix,]$Out,tabeic=tmpeic,npad=npad,stepmz=stepmz,mzdata=mzdata,chunk=chunk,verbose=FALSE,serial=FALSE)
  
}


#### wrapper for loadEICs
parseSampEIC<-function(matfile,tabeic,corrt=NULL,npad=3,stepmz=1/1000,verbose=TRUE,ncl=1,mzdata=FALSE,chunk=200){
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
  }
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  lsids=matfile$Sid[file.exists(matfile$In)]
  if(ncl==1){
    if(mzdata)   require("mzR")
    cat(" on 1 processor\n",sep="")
    re=list()
    for(ix in lsids){
      tmpeic=tabeic
      if("shppm"%in%names(matfile)) tmpeic$shppm=matfile[ix,]$shppm else tmpeic$shppm=0
      if("shrt"%in%names(matfile)) tmpeic$shrt=matfile[ix,]$shrt else tmpeic$shrt=0
      if(!is.null(corrt)) if(ix%in%rownames(corrt))
        tmpeic$shrt=tmpeic$shrt+approx(as.numeric(colnames(corrt)),corrt[ix,],tmpeic$rtap)$y
      
      re[[ix]]=parseOneSampEIC(matfile[ix,]$In,outfile=matfile[ix,]$Out,tabeic=tmpeic,npad=npad,stepmz=stepmz,mzdata=mzdata,chunk=chunk,verbose=verbose,serial=FALSE)
    }
    d1=proc.time()[3]
    cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(lsids),1)," secs per file\n",sep="")
    return(invisible(re))
  }
  cat(" on ",ncl," processors\n",sep="")
  sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile='sampeic.log')
  #if(sfIsRunning()) sfStop()
  if(mzdata) sfLibrary(mzR)
  sfLibrary(GRMeta)
  re=sfClusterApplyLB(lsids,.GRparseSampEICCl,matfile=matfile,corrt=corrt,tabeic=tabeic,npad=npad,stepmz=stepmz,mzdata=mzdata,chunk=chunk)
  sfStop()
  d1=proc.time()[3]
  cat(sapply(re,function(x) x[[1]]),sep="\n")
  lts=sapply(re,function(x) x[[2]])
  cat("Completed at ",date()," -> took ",round(d1-d0,1)," secs ",round((d1-d0)/length(lsids),1),"/",round(mean(lts),1)," secs per file\n",sep="")
  
  return(invisible(re))
  
}



parseOneSampEICold<-function(mzfi,tabeic,outfile=NULL,npad=3,stepmz=1/1000,verbose=TRUE,serial=TRUE){
  
  
  d0=proc.time()[3]
  if(!file.exists(mzfi)) return(NULL)
  if(!"shrt"%in%names(tabeic)) tabeic$shrt=0
  if(!"shppm"%in%names(tabeic)) tabeic$shppm=0
  ##########################################
  ## load data
  if(grepl('mzdata',mzfi)){
    aa=openMSfile(mzfi,backend="Ramp",verb=T)
    rts=header(aa)[,"retentionTime"]/60
    names(rts)=header(aa)[,1]
    pl <- mzR::peaks(aa)
    close(aa)
    m=do.call("rbind",lapply(1:length(pl),function(x) cbind(Scan=rep(x,nrow(pl[[x]])),pl[[x]])))
    rm("pl")
    m=m[order(m[,2]),]
    m=cbind(m,rts[m[,1]])
    colnames(m)=c("scan","mz","y","rt")
  }else load(mzfi)
  
  
  whichmz='mz';whichrt='rt'
  
  padrt=ifelse(npad<0,0,npad*median(diff(rts))*(npad+1))
  
  ##########################################
  ## filter the data file: wide mz
  rtmzrange=cbind(floor(tabeic[,'mzmin']*(1-tabeic[,'shppm']*10^-6-2*stepmz)/stepmz),
                  ceiling(tabeic[,'mzmax']*(1+tabeic[,'shppm']*10^-6+2*stepmz)/stepmz))
  
  ucode=apply(rtmzrange,1,function(x) x[1]:x[2])
  lcodecc=round(m[,whichmz]/stepmz)
  
  l2k=(lcodecc%in%unique(unlist(ucode)))
  mred=m[l2k,]
  lcodecc=lcodecc[l2k]
  if(verbose) cat("Reduction from ",nrow(m)," to ",nrow(mred)," ",sep="")
  rm(list="m")
  mred=mred[match(unlist(ucode),lcodecc),]
  mred=cbind(mred,eic=unlist(lapply(1:length(ucode),function(x) rep(x,length(ucode[[x]])))))
  if(npad<0) mred=mred[!is.na(mred[,1]),] 
  rm('ucode')
  ##############
  ## reduce mred
  mred=mred[which((mred[,whichrt]+tabeic[mred[,5],"shrt"]-tabeic[mred[,5],"rtmin"]-padrt)>=0),]
  mred=mred[which((mred[,whichrt]+tabeic[mred[,5],"shrt"]-tabeic[mred[,5],"rtmax"]+padrt)<=0),]
  l2rm=which((mred[,whichmz]*(1+tabeic[mred[,5],"shppm"]*10^-6)-tabeic[mred[,5],"mzmin"])<0)
  mred=mred[-l2rm,]
  l2rm=which((mred[,whichmz]*(1+tabeic[mred[,5],"shppm"]*10^-6)-tabeic[mred[,5],"mzmax"])>0)
  mred=mred[-l2rm,]
  mred<-mred[order(mred[,"eic"],mred[,"scan"]),]
  if(verbose) cat("to ",nrow(mred),sep="")
  #  cat("end at ",date(),sep="")
  
  ###################
  ## add missing scans
  if(npad>=0){
    lma=unlist(tapply(1:nrow(mred),mred[,"eic"],function(x)
      x[match(min(mred[x,'scan']):max(mred[x,'scan']),mred[x,'scan'])]))
    mred=mred[lma,]
    rm(list=c("lma"))
    for(i in which(is.na(mred[,"scan"]))){
      mred[i,"scan"]=mred[i-1,"scan"]+1
      mred[i,"eic"]=mred[i-1,"eic"]
    }
  }
  mred[which(is.na(mred[,'y'])),'y']=0
  mred[,'rt']=rts[as.character(mred[,'scan'])]
  if(verbose) cat(" plus ",nrow(mred),"\n",sep="")
  
  ########################
  ## add m/z and rt corr
  mred=cbind(mred,'rtcor'=mred[,'rt']+tabeic[mred[,'eic'],"shrt"])
  mred=cbind(mred,'mzcor'=mred[,'mz']*(1+tabeic[mred[,'eic'],"shppm"]*mred[,'mz']*10^-6))
  rownames(mred)=NULL
  
  d1=proc.time()[3]
  n=length(unique(mred[,'eic']))
  top=paste(mzfi," :"," -> ",round(d1-d0,1)," secs :",n,"/",nrow(tabeic)," (",round(100*n/nrow(tabeic),1),"%) EICs founds.",sep="")
  if(!is.null(outfile)) save(file=outfile,mred,rts,tabeic,mzfi)
  if(!serial) return(invisible(list(top,d1-d0)))
  if(verbose)   cat(top,sep="\n")
  return(invisible(mred))
  
}
