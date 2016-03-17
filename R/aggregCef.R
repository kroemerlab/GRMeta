aggregCEF<-function(lfiles,minrt=0,maxrt=Inf,maxmz=+Inf,minmz=1,ncl=1,verbose=TRUE){
  
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
    for(i in 1:length(lfiles)) acef[[i]]=.GRdoOneCef(lfiles[i],verbose=verbose)
  }
  if(ncl>1){
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile='logAggregCef')
    sfLibrary(GRMeta)
    #if(local) 
    #sfExport( ".GRdoOneCef", local=TRUE )
    acef=sfClusterApplyLB(lfiles,.GRdoOneCef,verbose=FALSE)
    sfStop()
  }  
  re=do.call("rbind",acef)
  re$samp=names(lfiles)[match(re$samp,lfiles)]
  re=re[which(re$rt>=minrt & re$rt<=maxrt),]
  re=re[which(re$mz>=minmz & re$mz<=maxmz),]
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(lfiles),1)," secs per file\n",sep="")
  return(invisible(re))
  
  
}

.GRdoOneCef<-function(cefi,verbose=TRUE){
  
  inOne<-function(ient,ix){
    
    lpks=xmlChildren(ient)$Location
    ## molecular feature location
    v=c(ix,paste(ix,".0",sep=""),(xmlAttrs(lpks)[c("rt","m","y","v")]),0,"M",0)
    names(v)=c("id","id2","rt","x","y","v","z","s","niso")
    
    ## individual peaks in the mol feat
    lpks=xmlChildren(xmlChildren(ient)$Spectrum)$MSPeaks
    mod=xmlAttrs(xmlChildren(xmlChildren(ient)$Spectrum)$MSDetails)
    dat2=t(xmlSApply(lpks,function(i) xmlAttrs(i)[c("rt","x","y","v","z","s")]))
    m <- gregexpr("\\+([1-9]+)$",dat2[,"s"])
    lma=regmatches(dat2[,"s"], m)
    niso=suppressWarnings(as.numeric(gsub("\\+","",lma)))
    niso[is.na(niso)]=0
    id2=paste(ix,1:length(niso),sep=".")
    dat2=suppressWarnings(data.frame(id=ix,id2=id2,dat2,niso=niso,stringsAsFactors=F,row.names=NULL))
    rownames(dat2)=NULL
    dat2=data.frame(samp=NA,rbind(v,dat2),stringsAsFactors=FALSE)
    for(i in c("rt","x","y","v","z","niso")) dat2[,i]=as.numeric(gsub(",",".",dat2[,i]))
    names(dat2)=c("samp", "mfid","mfid2","rt","mz","Area","Height","xM","ip","iso")
    return(dat2)
  }
  
  if(verbose) cat(cefi,":")
  doc <- xmlParse(cefi, useInternal=TRUE)
  xml=xmlChildren(xmlChildren(doc)[[1]])[[1]]
  adat=lapply(1:xmlSize(xml),function(i) suppressWarnings(inOne(xml[[i]],i)))
  cefdat=do.call("rbind",adat)
  rm(list=c("adat","xml"))
  free(doc)
  cefdat$samp=cefi
  if(verbose) cat(" ",length(unique(cefdat$mfid))," mol feat./",nrow(cefdat)," entries\n",sep="")
  return(cefdat)
  
}
