
parseOneSampEIC<-function(mzfi,tabeic,outfile=NULL,npad=3,stepmz=1/1000,mzdata=FALSE,addMZ=1,addRT=0.1,extraRT=NA,extraPPM=NA,keepM=FALSE,verbose=TRUE,serial=TRUE,chunk=200){
  
  
  d0=proc.time()[3]
  if(!file.exists(mzfi)) return(NULL)
  # if(!"shrt"%in%names(tabeic)) tabeic$shrt=0
  # if(!"shppm"%in%names(tabeic)) tabeic$shppm=0
  ##########################################
  ## load data
  if(mzdata){
    aa=openMSfile(mzfi,backend="Ramp",verb=T)
    rts=mzR::header(aa)[,"retentionTime"]/60
    names(rts)=mzR::header(aa)[,1]
    pl <- mzR::peaks(aa)
    mzR::close(aa)
    m=do.call("rbind",lapply(1:length(pl),function(x) cbind(Scan=rep(x,nrow(pl[[x]])),pl[[x]])))
    rm("pl")
    m=m[order(m[,2]),]
    m=cbind(m,rts[m[,1]])
    colnames(m)=c("scan","mz","y","rt")
  }else load(mzfi)
  
 
  whichmz='mz'
  if(!is.na(addMZ)){
    if(!"mzcor"%in%colnames(m)) m=cbind(m,mzcor=m[,'mz']*(1+addMZ*10^-6)) else m[,'mzcor']=m[,'mz']*(1+addMZ*10^-6)
    whichmz='mzcor'
    m[,'mzcor']=round(m[,'mzcor'],7)
  }
  
  whichrt='rt'
  if(!is.na(addRT)){
    if(!"rtcor"%in%colnames(m)) m=cbind(m,rtcor=m[,'rt']+addRT) else m[,'rtcor']=m[,'rt']+addRT
    whichrt='rtcor'
    m[,'rtcor']=round(m[,'rtcor'],ceiling(abs(log10(median(diff(rts))))+3))
  }
  mnames=colnames(m)
  
  tabeic=tabeic[order(tabeic[,'mzmin']),]
  
  
#  padrt=ifelse(npad<0,0,npad*median(diff(rts))*(npad+1))
  
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
    mzmin=itab[,'mzmin']
    mzmax=itab[,'mzmax']
    if(!is.na(extraPPM)){
      mzmin=mzmin-itab[,"mzap"]*extraPPM*10^-6
      mzmax=mzmax+itab[,"mzap"]*extraPPM*10^-6
    }
    rtmin=itab[,'rtmin']
    rtmax=itab[,'rtmax']
#    print(rtmin)
    if(!is.na(extraRT)){
      rtmin=rtmin-extraRT
      rtmax=rtmax+extraRT
    }
    names(rtmin)=names(rtmax)=names(mzmin)=names(mzmax)=itab[,"Id"]
#    print(rtmin)
    rtmzrange=cbind(floor((mzmin-3*stepmz)/stepmz),
                    ceiling((mzmax+3*stepmz)/stepmz))
    l=which(lcode>=min(rtmzrange) & lcode<=max(rtmzrange))
    if(length(l)==0) next
    v=lcode[l]
    mred=m[l,,drop=F]
 #   print(dim(mred))
#    system.time(rangmz<-apply(rtmzrange,1,function(x) xcms:::findRange(v,x[1],x[2])))
    #ucode=apply(t(rangmz),1,function(x) x[1]:x[2])
    system.time(ucode<-apply(rtmzrange,1,function(x) which(v>=x[1] & v <=x[2])))
     mred=mred[unlist(ucode),,drop=F]
 #    print(c(3,dim(mred)))
     if(nrow(mred)>0){
     mred=cbind(mred,eic=unlist(lapply(1:length(ucode),function(x) rep(x,length(ucode[[x]])))))
#     print(dim(mred))
     mred=mred[which((mred[,whichrt]-rtmin[mred[,'eic']])>=0),,drop=FALSE]
    mred=mred[which((mred[,whichrt]-rtmax[mred[,'eic']])<=0),,drop=FALSE]
#    print(dim(mred))
    mred=mred[which((mred[,whichmz]-mzmin[mred[,'eic']])>=0),,drop=FALSE]
    mred=mred[which((mred[,whichmz]-mzmax[mred[,'eic']])<=0),,drop=FALSE]
    mred=cbind(mred,ineic=0)
    linrt=(mred[,whichrt]-itab[mred[,'eic'],"rtmin"])>=0 & (mred[,whichrt]-itab[mred[,'eic'],"rtmax"])<=0
    linmz=(mred[,whichmz]-itab[mred[,'eic'],"mzmin"])>=0 & (mred[,whichmz]-itab[mred[,'eic'],"mzmax"])<=0
    mred[,"ineic"]=(linrt & linmz)*1
    mred[,"eic"]=lins[mred[,"eic"]]
     } else{mred=matrix(nrow=0,ncol=length(mnames)+1);colnames(mred)=c(mnames,'eic')}
#     print(c(5,dim(mred)))
     amred[[ilins]]=mred
    rm(list="mred")
  }
  rm(list=c('lcode','llcode'))
  if(!keepM) rm(list="m")
  mred=do.call("rbind",amred)
  rm(list=c('amred'))
  mred<-mred[order(mred[,"eic"],mred[,"scan"]),]
  mred[which(is.na(mred[,'y'])),'y']=0
 # mred[,'rt']=rts[as.character(mred[,'scan'])]
  if(verbose) cat(" plus ",nrow(mred),"\n",sep="")
  
  rownames(mred)=NULL  
  attr(mred,"eic")<-tabeic
  attr(mred,"file")<-mzfi
  attr(mred,"rts")<-rts
  attr(mred,"addRT")<-addRT
  attr(mred,"addMZ")<-addMZ

  d1=proc.time()[3]
  n=table(mred[,'eic'])
  top=paste(mzfi," :"," -> ",round(d1-d0,1)," secs :",length(n),"/",nrow(tabeic),": ",round(mean(n),3)," scan "," (",round(100*length(n)/nrow(tabeic),1),"%) EICs founds.",sep="")
  if(!is.null(outfile))
    if(keepM) save(file=outfile,m,mred) else save(file=outfile,mred)
  if(!serial) return(invisible(list(top,d1-d0)))
  if(verbose)   cat(top,sep="\n")
  return(invisible(mred))
  
}


.GRparseSampEICCl<-function(ix,matfile,tabeic,stepmz,mzdata,keepM,extraRT,extraPPM,chunk,addCor=TRUE){
  addMZ=ifelse("shppm"%in%names(matfile),matfile[ix,]$shppm,ifelse(addCor,0,NA))
  addRT=ifelse("shrt"%in%names(matfile),matfile[ix,]$shrt,ifelse(addCor,0,NA))
  
 parseOneSampEIC(matfile[ix,]$In,tabeic=tabeic,outfile=matfile[ix,]$Out,stepmz=stepmz,
                mzdata=mzdata,addRT=addRT,addMZ=addMZ,keepM=keepM,extraRT=extraRT,extraPPM=extraPPM,
                chunk=chunk,verbose=FALSE,serial=FALSE)
# parseOneSampEIC<-function(mzfi,tabeic,outfile=NULL,npad=3,stepmz=1/1000,mzdata=FALSE,addMZ=NA,addRT=NA,verbose=TRUE,serial=TRUE,chunk=200){
   
}


#### wrapper for loadEICs
parseSampEIC<-function(matfile,tabeic,stepmz=1/1000,keepM=TRUE,extraRT=NA,extraPPM=NA,verbose=TRUE,ncl=1,mzdata=FALSE,addCor=TRUE,chunk=200){
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,nrow(matfile),parallel:::detectCores()))
  }
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  lsids=matfile$Sid[file.exists(matfile$In)]
  if(ncl==1){
    if(mzdata)   require("mzR")
    cat(" on 1 processor\n",sep="")
    re=list()
    for(ix in lsids){
      addMZ=ifelse("shppm"%in%names(matfile),matfile[ix,]$shppm,ifelse(addCor,0,NA))
      addRT=ifelse("shrt"%in%names(matfile),matfile[ix,]$shrt,ifelse(addCor,0,NA))
      
      re[[ix]]=parseOneSampEIC(matfile[ix,]$In,outfile=matfile[ix,]$Out,tabeic=tabeic,
            stepmz=stepmz,mzdata=mzdata,addRT=addRT,addMZ=addMZ,keepM=keepM,extraRT=extraRT,extraPPM=extraPPM,chunk=chunk,verbose=verbose,serial=TRUE)
    }
    d1=proc.time()[3]
    cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(lsids),1)," secs per file\n",sep="")
    return(invisible(re))
  }
  cat(" on ",ncl," processors\n",sep="")
  sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile=paste('SampEic',format(Sys.time(), "%y-%d-%b-%H:%M"),'.log',sep=""))
  #if(sfIsRunning()) sfStop()
  if(mzdata) sfLibrary(mzR)
  sfLibrary(GRMeta)
  re=sfClusterApplyLB(lsids,.GRparseSampEICCl,matfile=matfile,tabeic=tabeic,stepmz=stepmz,mzdata=mzdata,keepM=keepM,
                     extraRT=extraRT,extraPPM=extraPPM,chunk=chunk,addCor=addCor)
#  .GRparseSampEICCl<-function(ix,matfile,tabeic,npad,stepmz,mzdata,chunk,addCor=TRUE){
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
    rts=mzR::header(aa)[,"retentionTime"]/60
    names(rts)=mzR::header(aa)[,1]
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
      x[match(min(mred[x,'scan']):max(mred[x,'scan']),mred[x,'scan'])])) ### dodgy implementation: pb if several mz for same scan??
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
