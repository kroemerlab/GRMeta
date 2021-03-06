##############
## 

matchmzSet<-function(obj,Analyte=NULL,annotdb=NULL,ipdb=NULL,lIP=NULL,
                   dppm=20,mz2use="MZ",mass2use="Mass",
                   annotinf2add=NULL,objinf2add=NULL,groupby=NULL,chunk=200,ncl=1,collsep=";;"){
  
  if(inherits(obj, "metaboSet")){
    cat("Matching from a metaboSet\n")
    dfIn=obj$Annot
    if(is.null(Analyte)) Analyte=obj$Analyte
    lanalytes=Analyte[Analyte%in%obj$Analyte]
  }else if(inherits(obj,"data.frame")){
    cat("Matching from a data.frame\n")
    dfIn=obj
    if(is.null(rownames(dfIn))) rownames(dfIn)=paste("Ana",1:nrow(dfIn),sep="")
    if(any(duplicated(rownames(dfIn)))) rownames(dfIn)=paste("Ana",1:nrow(dfIn),sep="")
    if(is.null(Analyte)) Analyte=rownames(dfIn)
    lanalytes=Analyte[Analyte%in%rownames(dfIn)]
  #  print(names(dfIn))
  }else if(inherits(obj,"numeric") & is.null(dim(obj))){
    cat("Matching from a vector\n")
    if(is.null(names(obj))) names(obj)=paste("Ana",1:length(obj),sep="")
    dfIn=data.frame(MZ=obj,stringsAsFactors=FALSE)
    rownames(dfIn)=names(obj)
    if(is.null(Analyte)) Analyte=rownames(dfIn)
    lanalytes=Analyte[Analyte%in%rownames(dfIn)]
  }else stop("Provide metaboSet, data.frame or a vector")
  
  if(is.null(ipdb)){
    data(IPDB)
    ipdb=IPDB
#    rm('IPDB')
  }
  if(is.null(annotdb)){
    data(AnnotationDB)
    annotdb=AnnotationDB
#    rm('AnnotationDB')
  }
  
  
  if(length(lanalytes)==0) stop("No analyte found in the object\n")
  ######################
  if(!mz2use%in%names(dfIn)) stop(paste("m/z variable",mz2use,"not found in the object\n"))
  cmz=as.numeric(dfIn[lanalytes,mz2use])
  names(cmz)=lanalytes
  cmz=cmz[!is.na(cmz)]
  if(length(cmz)==0) stop(paste("nothing to match in the m/z variable",mz2use,"\n"))

  ######################
  if(!mass2use%in%names(annotdb)) stop(paste("Mass variable",mass2use,"not found in the annotation dataframe\n"))
  lentries=which(!is.na(as.numeric(annotdb[,mass2use])))
  if(length(lentries)<1) stop("nothing to match in the annotation dataframe\n")
  annotdb=annotdb[lentries,,drop=F]
  lmz=annotdb[,mass2use]
  
  
  ladd=unique(which(ipdb$Id%in%lIP | ipdb$Name%in%lIP | ipdb$Set%in%lIP))
  if(length(ladd)==0)
    stop(paste("IPs invalid - must be in :\n",
               "  *positive mode: ",paste(IPDB$Name[IPDB$Charge>0],collapse=" "),"\n",
               "  *negative mode: ",paste(IPDB$Name[IPDB$Charge<0],collapse=" "),"\n",sep=""))
  cat("Matching on: ",ipdb$Name[ladd],"\n",sep=" ")
  mmz=matrix(sapply(ladd,function(i) (ipdb[i,]$xM*lmz+ipdb[i,]$adMass)/abs(ipdb[i,]$Charge)),nrow=length(lmz))
  colnames(mmz)=ipdb$Name[ladd]
  ######################
  if(length(cmz)<=chunk) lx=list(names(cmz)) else{
  ngrps=max(1,round(length(cmz)/chunk))
  lst=seq(1,length(cmz),length.out=ngrps+1)
  lx=suppressWarnings(split(names(cmz),1:ngrps))
  }
  ### parallel? lx/mmz/cmz/dppm
  
  .infctMZM<-function(cx,mmz,cmz,dppm){
    
    allm=list()
    for(i in 1:ncol(mmz)){
      mdppm=sweep(outer(cmz[cx],mmz[,i],"-"),1,cmz[cx],"/")*10^6
      lmatch=which(abs(mdppm)<=dppm,arr=T)
      if(length(lmatch)>0) 
        allm[[i]]=data.frame(Analyte=cx[lmatch[,1]],IP=i,DPPM=round(mdppm[lmatch],3),Entry=lmatch[,2],stringsAsFactors=F)
    }
    if(length(allm)>0) return(do.call("rbind",allm))
    return(NULL)
    
  }
  

  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  
  if(ncl==1 | length(lx)==1){
    cat(" on 1 processor\n",sep="")
    allr=list()
    for(k in 1:length(lx)){
      cx=lx[[k]]
      cat(ifelse(k%%10==0,"X","."))
      allr[[k]]=.infctMZM(cx,mmz,cmz,dppm)
    }
  }
  if(ncl>1 & length(lx)>1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile='loggetTIC')
    sfExport( ".infctMZM", local=TRUE )
    allr=sfClusterApplyLB(lx,.infctMZM,mmz,cmz,dppm)
    sfStop()
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs \n",sep="")
  
  allr=allr[which(!sapply(allr,is.null))]
  cat("\n")
  if(length(allr)==0) return(NULL)
  matchres=do.call("rbind",allr)

  ############
  matchres$IP=colnames(mmz)[matchres$IP]
  rownames(matchres)=NULL
  if(!is.null(annotinf2add)){
    annotinf2add=unique(c(annotinf2add[annotinf2add%in%names(annotdb)],groupby[groupby%in%names(annotdb)]))
    if(length(annotinf2add)>0){
      cat("* adding from annotdb: ",annotinf2add,"\n",sep=" ")
      for(i in annotinf2add) matchres[,i]=annotdb[matchres$Entry,i]
    }
  }
  if(!is.null(objinf2add)){
    objinf2add=unique(objinf2add[objinf2add%in%names(dfIn) & !objinf2add%in%names(matchres)])
 #   print(objinf2add)
    if(length(objinf2add)>0){
      cat("* adding from obj: ",objinf2add,"\n",sep=" ")
      for(i in objinf2add) matchres[,i]=dfIn[matchres$Analyte,i]
    }
  }
  
  matchres$Entry=lentries[matchres$Entry]
  ############
  iDPPM=matchres$DPPM
  if(!is.null(groupby)){
    groupby=groupby[groupby%in%names(matchres)]
    if(length(groupby)>0){
      cat("* grouping by: ",groupby,"\n",sep=" ")
      grpfac=apply(matchres[,c("Analyte","IP",groupby)],1,function(x) paste(as.character(x),collapse=".;."))
      lnam=names(matchres)
      iDPPM=round(tapply(iDPPM,grpfac,median),3)
      matchres=lapply(names(matchres),function(x) tapply(matchres[,x],grpfac,function(y) paste(unique(y),collapse=collsep)))
      matchres=data.frame(do.call("cbind",matchres),stringsAsFactors=FALSE)
      names(matchres)=lnam
      if(any(grepl(collsep,matchres$DPPM))) matchres$DPPMmed=iDPPM
    }
  }
  ############  
  matchres=matchres[order(factor(matchres$Analyte,levels=lanalytes),factor(matchres$IP,levels=colnames(mmz)),abs(iDPPM)),]
  rownames(matchres)=NULL
  
  cat("  -> ",length(unique(matchres$Analyte))," unique analytes\n",sep="")
  cat("  -> ",length(unique(matchres$Entry))," unique entries\n",sep="")
  summ=c(median(matchres$DPPM),mad(matchres$DPPM),mean(matchres$DPPM),sd(matchres$DPPM),range(matchres$DPPM))
  cat("  -> ",sprintf("Med= %.2f(%.2f) Mean=%.2f(%.2f) [%.2f;%.2f]",summ[1],summ[2],summ[3],summ[4],summ[5],summ[6])," unique entries\n",sep="")
  
  return(matchres)
}
