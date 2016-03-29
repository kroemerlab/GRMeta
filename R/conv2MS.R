.inGRconvMS<-function(lfiles,lsids,
                      ldatas=c(MZ="MZap",Area="Area.scan",AreaSm="Area.scansm",
                               Height="HE.ap",RT="RT.ap",RTsm="RT.apsm",Coda="Coda",NCons="NCons")){
  alleics=alldata=list()
  for(ifi in 1:length(lfiles)){
    system.time(load(lfiles[ifi]))
    rm("dfeic")
    alldat=lapply(alldat,function(x) x[match(lsids,rownames(x)),,drop=F])
    alldat=alldat[names(alldat)%in%ldatas]
    for(i in names(alldat)) if(is.null(alldata[[i]])) alldata[[i]]=alldat[[i]] else alldata[[i]]=cbind(alldata[[i]],alldat[[i]])
    alleics[[ifi]]=eicstats
  }
  alldata=alldata[!sapply(alldata,is.null)]
  eicsdat=do.call("rbind",alleics)
  rm(list=c("alleics",'alldat','eicstats'))
  
  leics=eicsdat$PkId
  alldata=lapply(alldata,function(x) x[,match(leics,colnames(x)),drop=F])
  names(alldata)=names(ldatas[match(names(alldata),ldatas)])
  for(k in names(alldata)) dimnames(alldata[[k]])=list(lsids,leics)
  return(list(dat=alldata,eic=eicsdat))
}

conv2metaboSet<-function(lfiles,Meta,File,method="prof",
                      ldatas=c(MZ="MZap",Area="Area.scan",AreaSm="Area.scansm",
                               Height="HE.ap",RT="RT.ap",RTsm="RT.apsm",Coda="Coda",NCons="NCons"),chunk=10,ncl=1){
  
  if(any(!c("Sid","sType")%in%names(Meta))) stop('Meta must contain: Sid/sType')
  if(any(!c("Sid","File","Date")%in%names(File))) stop('File must contain: Sid/File/Date')
  
  rownames(Meta)=Meta$Sid
  lsids=Meta$Sid
  
  File=File[match(lsids,File$Sid),]
  rownames(File)=File$Sid=lsids
  
  lfiles=lfiles[file.exists(lfiles)]
  if(length(lfiles)==0) stop('Files not found')
  
  if(length(lfiles)<=chunk) llx=list(lfiles) else  llx=split(lfiles, ceiling(seq_along(lfiles)/chunk))
  
  
  ##########################
  d0=proc.time()[3]
  cat("Processing ",length(lfiles)," files: starting at ",date(),sep="")
  
  if(ncl==1 | length(llx)==1){
    cat(" on 1 processor\n",sep="")
    allr=list()
    for(k in 1:length(llx)){
      cx=lx[[k]]
      cat(ifelse(k%%10==0,"X","."))
      system.time(allr[[k]]<-.inGRconvMS(llx[[k]],lsids,ldatas))
    }
  }
  if(ncl>1 & length(llx)>1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfExport( ".inGRconvMS", local=TRUE )
    allr=sfClusterApplyLB(llx,.inGRconvMS,lsids,ldatas)
    sfStop()
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs \n",sep="")
  
  ##########################
  alldata=list()
  for(x in allr)
    for(i in names(x[[1]])) if(is.null(alldata[[i]])) alldata[[i]]=x[[1]][[i]] else alldata[[i]]=cbind(alldata[[i]],x[[1]][[i]])
  alleics=list()
  for(x in 1:length(allr)) alleics[[x]]=allr[[x]]$eic
  
  eicsdat=do.call("rbind",alleics)

  ########## prep annot
  imz=apply(alldata$MZ,2,median,na.rm=T)
  irt=apply(alldata$RT,2,median,na.rm=T)
  newn=sprintf("%.4f@%.2f-%s",imz,irt,method)
  ldups=names(which(table(newn)>1))
  newn[newn%in%ldups]=sprintf("%.5f@%.3f-%s",imz,irt,method)[newn%in%ldups]
  ldups=names(which(table(newn)>1))
  for(i in ldups){
    l=which(newn==i)
    newn[l]=sprintf("%.5f-D%d@%.3f-%s",imz[l],1:length(l),irt[l],method)
  }
  Annot=data.frame(Analyte=newn,MetName=NA,IsSTD=FALSE,OriginalName=eicsdat$PkId,LevelAnnot=4,Method=method,MZ=round(imz,6),RT=round(irt,4),stringsAsFactors=F)
  rownames(Annot)=newn
  for(i in names(alldata)) colnames(alldata[[i]])=newn
  
  #########
  obj=list(Method=method,Sid=lsids,Analyte=newn,Annot=Annot,Meta=Meta,File=File,Eic=list(Path=getwd(),File=eicsdat),Data=alldata)
  class(obj)=append(class(obj),"metaboSet")
print(obj)
  invisible(obj)
}
