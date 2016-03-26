.GRgetTIC<-function(mzfi,verbose=TRUE){
  require(mzR)
  if(!file.exists(mzfi)){
    cat(paste(mzfi,"does not exist"))
    return(list(matrix(c(1,NA,NA),ncol=3),mzfi))
  }
  aa=openMSfile(mzfi,backend="Ramp",verb=T)
  rts=header(aa)[,"retentionTime"]/60
  names(rts)=header(aa)[,1]
  pl <- mzR::peaks(aa)
  mzR::close(aa)
  m=do.call("rbind",lapply(1:length(pl),function(x) cbind(Scan=rep(x,nrow(pl[[x]])),pl[[x]])))
  rm("pl")
  if(verbose) cat("* ",mzfi,": mz=",round(min(m[,2]),5),"-",round(max(m[,2]),5),
                  ": rt=",round(min(rts),2),"-",round(max(rts),2),"\n")
  tic=tapply(m[,3],m[,1],sum,na.rm=T)
  rm("m")
  list(cbind(scan=as.numeric(names(tic)),rt=rts[as.numeric(names(tic))],tic=tic),mzfi)
}


.GRmgetTIC<-function(lfiles,ncl=1,verbose=(ncl==1),mergeTIC=TRUE){
  require(mzR)
  
  if(is.null(names(lfiles))) names(lfiles)=1:length(lfiles)
  lfiles=lfiles[file.exists(lfiles)]
  if(length(lfiles)==0) return(NULL)
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
  }
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  
  if(ncl==1){
    acef=list()
    cat(" on 1 processor\n",sep="")
    for(i in 1:length(lfiles)) acef[[i]]=.GRgetTIC(lfiles[i],verbose=verbose)
  }
  if(ncl>1){
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile='loggetTIC')
    sfLibrary(GRMeta)
    sfLibrary(mzR)
    #if(local) 
    #sfExport( ".GRdoOneCef", local=TRUE )
    acef=sfClusterApplyLB(lfiles,.GRgetTIC,verbose=FALSE)
    sfStop()
  }
  names(acef)=names(unlist(sapply(acef,function(x) x[2])))
  re=lapply(acef,function(x) x[[1]])
  if(mergeTIC){
  rscan=range(sapply(re,function(x) range(x[,1],na.rm=T)),na.rm=T)
  lsc=min(rscan):max(rscan)
  re=lapply(re,function(x) x[match(lsc,x[,1]),,drop=F])
  rt=sapply(re,function(x) x[,2,drop=F])
  inty=sapply(re,function(x) x[,3,drop=F])
  rownames(rt)=rownames(inty)=lsc
  re=list(rt=rt,inty=inty)
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(lfiles),1)," secs per file\n",sep="")
  return(invisible(re))

  }


