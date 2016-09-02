aggregCEF<-function(lfiles,minrt=0,maxrt=Inf,maxmz=+Inf,minmz=1,ncl=1,type="MS",verbose=TRUE){
  
  if(!type%in%c('MS','AIMS')) stop('Wrong type')
  
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
    for(i in 1:length(lfiles)) acef[[i]]=.GRdoOneCef(lfiles[i],type=type,verbose=verbose)
  }
  if(ncl>1){
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile='logAggregCef')
    sfLibrary(GRMeta)
    #if(local) 
    #sfExport( ".GRdoOneCef", local=TRUE )
    acef=sfClusterApplyLB(lfiles,.GRdoOneCef,type=type,verbose=FALSE)
    sfStop()
  }  
    re=do.call("rbind",acef)
    re$samp=names(lfiles)[match(re$samp,lfiles)]
   if(type=="MS"){
     re=re[which(re$rt>=minrt & re$rt<=maxrt),]
    re=re[which(re$mz>=minmz & re$mz<=maxmz),]
   }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(lfiles),1)," secs per file\n",sep="")
  return(invisible(re))
  
  
}

.GRdoOneCef<-function(cefi,type="MS",verbose=TRUE){
  if(!type%in%c('MS','AIMS')) stop('Wrong type')
  inOneAION<-function(ient,ix){
    
    ient2=xmlChildren(ient)
    
    lmspk=which(names(ient2)=="Spectrum")
    infos=pks=list()
    for(j in lmspk){
      k=j-lmspk[1]+1
      x=xmlChildren(ient2[[j]])
      infos[[k]]=(c(k,xmlAttrs(x$MSDetails)[c("p","ce","fv")],
                    as.numeric(xmlAttrs(x$RTRanges[[1]])),
                    xmlValue(x$MzOfInterest)))
      pks[[k]]=cbind(t(xmlSApply(x$MSPeaks,function(x) as.numeric(xmlAttrs(x)[c("x","y")]))),k)
    }
    infos=data.frame(do.call("rbind",infos),stringsAsFactors=FALSE)
    names(infos)=c("id","mod","ce","fv","rtmin","rtmax","parent")
    
    #######
    infos$ce=as.numeric(gsub("[A-Za-z]","",infos$ce))
    infos$ce[is.na(infos$ce)]=0
    infos$fv=as.numeric(gsub("[A-Za-z]","",infos$fv))
    infos$fv[is.na(infos$fv)]=0
    infos$rtmin=as.numeric(infos$rtmin)
    infos$rtmax=as.numeric(infos$rtmax)
    infos$parent=as.numeric(infos$parent)
    
    #######
    nid=sprintf("%.4f@%.3f",infos$parent[1],(infos$rtmin[1]+infos$rtmax[1])/2)
    
    #######
    pks=do.call("rbind",pks)
    colnames(pks)=c("mz","y","id2")
    rownames(pks)=1:nrow(pks)
    pks=data.frame(samp=NA,mfid=nid,mfid2=paste(ix,pks[,3],sep="."),pks[,1:2],infos[pks[,3],c("ce","parent","rtmin","rtmax")],stringsAsFactors=FALSE)
    return(Pks=pks)
    #return(list(Id=nid,Infos=infos,Pks=pks))
  }
  
  inOne<-function(ient,ix){
    
    lpks=xmlChildren(ient)$Location
    ## molecular feature location
    v=c(ix,paste(ix,".0",sep=""),(xmlAttrs(lpks)[c("rt","m","y","v")]),0,"M",0)
    names(v)=c("id","id2","rt","x","y","v","z","s","niso")
    
    ## individual peaks in the mol feat
    lpks=xmlChildren(xmlChildren(ient)$Spectrum)$MSPeaks
    mod=xmlAttrs(xmlChildren(xmlChildren(ient)$Spectrum)$MSDetails)
    dat2=t(xmlSApply(lpks,function(i) xmlAttrs(i)[c("rt","x","y","v","z","s")]))
    m <- gregexpr("\\+([0-9]+)$",dat2[,"s"])
    lma=regmatches(dat2[,"s"], m)
    niso=suppressWarnings(as.numeric(gsub("\\+","",lma)))
    niso[is.na(niso)]=0
    id2=paste(ix,1:length(niso),sep=".")
    dat2=suppressWarnings(data.frame(id=ix,id2=id2,dat2,niso=niso,stringsAsFactors=F,row.names=NULL))
    rownames(dat2)=NULL
    dat2=data.frame(samp=NA,rbind(v,dat2),stringsAsFactors=FALSE)
    for(i in c("rt","x","y","v","z","niso")) dat2[,i]=as.numeric(gsub(",",".",dat2[,i]))
    names(dat2)=c("samp", "mfid","mfid2","rt","mz","Height","Area","z","ip","iso")
    dat2$xM=1
    dat2$xM[which(dat2$ip=="M")]=0
    dat2$xM[grep("^[0-9]+M",dat2$ip)]=as.numeric(gsub("^([0-9]+)M.*","\\1",dat2$ip[grep("^[0-9]+M",dat2$ip)]))
    return(dat2)
  }
  
  if(verbose) cat(cefi,":")
  doc <- xmlParse(cefi, useInternal=TRUE)
  xml=xmlChildren(xmlChildren(doc)[[1]])[[1]]
  if(type=="MS") adat=lapply(1:xmlSize(xml),function(i) suppressWarnings(inOne(xml[[i]],i)))
  if(type=="AIMS") adat=lapply(1:xmlSize(xml),function(i) suppressWarnings(inOneAION(xml[[i]],i)))
  cefdat=do.call("rbind",adat)
  rm(list=c("adat","xml"))
  free(doc)
  cefdat$samp=cefi
  if(verbose) cat(" ",length(unique(cefdat$mfid))," mol feat./",nrow(cefdat)," entries\n",sep="")
  return(cefdat)
  
}
