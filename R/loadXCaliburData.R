loadXCaliburData<-function(ifile,ofile=NULL,params=list()){
  
  
  paramsvals <- paramsParsing()
  if (!missing(params)) paramsvals[names(params)] <- params
  params=paramsvals
  
  .infct<-function(rows){
    cells <- getCells(rows)
    res <- lapply(cells, getCellValue)
    ishidx=do.call("rbind",lapply(names(res),function(i) c(strsplit(i,"\\.")[[1]],res[[i]])))
    nr=max(as.numeric(ishidx[,1]))
    nc=1:max(as.numeric(ishidx[,2]))
    shmat=do.call("cbind",tapply(1:nrow(ishidx),factor(ishidx[,2],nc),function(x) ishidx[x[match(1:nr,ishidx[x,1])],3]))
    metnam=shmat[3,1]
    shmat=shmat[which(shmat[,1]=="Filename"):nrow(shmat),,drop=F]
    colnames(shmat)=shmat[1,]
    shmat=shmat[-1,]
    
    whichfile=which(colnames(shmat)=='Filename')
    fnames=shmat[,whichfile]
    
    
    l2use=which(!is.na(fnames) & !fnames%in%c("Created By:","User Name","Thermo"))
    whichdir=which(colnames(shmat)=='Exp Method')
    if(length(whichdir)) fnames=paste0(shmat[,whichdir],"\\",fnames)
    
    whichdate=which(colnames(shmat)=='Acq Date')
    dts=shmat[l2use,whichdate]
    ldats=c("Area","Height",  "RT","Del RT", "S/N","Start Time","End Time","Start Height","End Height")
    idats=suppressWarnings(apply(shmat[l2use,colnames(shmat)%in%ldats],2,as.numeric))
    rownames(idats)=names(dts)=fnames[l2use]
    list(idats,dts,metnam)
  }
  
#  params= paramsParsing()
  
  wb     <- loadWorkbook(ifile)
  sheets <- getSheets(wb)
  sheets=sheets[which(!names(sheets)%in%c( "Component", "mdlCalcs" ))]
  if(length(sheets)==0) stop('No metabolites found???')
  
  adats=lapply(names(sheets),function(isheet) .infct(getRows(sheets[[isheet]])))
  
  #################
  lumetnams=c()
  lfi=unique(unlist(lapply(adats,function(x) names(x[[2]]))))
  ldatatype=unique(unlist(lapply(adats,function(x) colnames(x[[1]]))))
  for(i in 1:length(adats)){
    lumetnams=c(lumetnams,adats[[i]][[3]])
    ldatatype=unique(c(ldatatype,colnames(adats[[i]][[1]])))
    adats[[i]][[1]]=adats[[i]][[1]][match(lfi,rownames(adats[[i]][[1]])),match(ldatatype,colnames(adats[[i]][[1]])),drop=F]
    adats[[i]][[2]]=as.numeric(adats[[i]][[2]][lfi])
  }
  #########################################################################################################################
  
  #########################################################################################################################
  dts=apply(do.call("cbind",lapply(adats,function(x) x[[2]])),1,mean,na.rm=T)
  nams=sapply(strsplit(lfi,"\\\\"),function(x) rev(x)[1])
  
  # 
  #   if(!is.null(params$NameCol)){
  #     whichname= which(as.matrix(tab[1,])[1,]==params$NameCol )
  #     if(length(whichname)==0) params$NameCol=NULL else  nams=(as.matrix(tab[,whichname]))[-1,1]
  #   } 
  #   if(is.null(params$NameCol))  nams=gsub("\\.[dD]$","",filenam)
  if(!is.null(params$NameClean)) for(i in params$NameClean) nams=gsub(i,"",nams)
  styp=rep(NA,length(nams))
  if(length(grep(params$regTypes,nams))>0)
    styp[grep(params$regTypes,nams)]=gsub(params$regTypes,"\\1",nams)[grep(params$regTypes,nams)]
  
  if(any(is.na(styp))){
    cat(sum(is.na(styp)),"unknown sample type:\n", nams[which(is.na(styp))],"\n",sep=" ")
    styp[is.na(styp)]="Unknown"
  }
  sid=nams
  while(any(table(sid)>1)){
    ldups=names(which(table(sid)>1))
    cat(length(ldups),"duplicated sample names:\n", ldups,"\n",sep=" ")
    for(i in ldups) sid[which(sid==i)[-1]]=paste(sid[which(sid==i)[-1]],letters[1:(sum(sid==i)-1)],sep="")
  }
  
  metainfos=data.frame(Sid=sid,sType=styp,InjOrder=order(order(dts)),stringsAsFactors=FALSE)
  rownames(metainfos)=metainfos$Sid
  fileinfos=data.frame(File=lfi,Date=chron(dts),Name=nams,Sid=sid,stringsAsFactors=FALSE)
  rownames(fileinfos)=fileinfos$Sid
  if(!is.null(params$Batch)) fileinfos$Batch=params$Batch
  
  
  #########################################################################################################################
  allmat=lapply(ldatatype,function(x){
    xm=do.call("cbind",lapply(adats,function(y){
      re=y[[1]][,x,drop=F]
      colnames(re)=y[[3]]
      re
    }))
    rownames(xm)=metainfos$Sid
    xm
  })
  names(allmat)=ldatatype
  
  #########
  rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
  annot=data.frame(Analyte=NA,MetName=lumetnams,IsSTD=FALSE,IsISO=FALSE,RT=rtmed,LevelAnnot=1,stringsAsFactors = F)
  newnam=oldnam=lumetnams
  if(params$checkNams){
    lnnunk=which(!grepl("_UNK$",newnam))
    newnam[lnnunk]=cleanMetaboNames(oldnam[lnnunk],RegExpr = NA,Syno = NA)$newnam
    #   print(newnam)
    annot$MetName=newnam
    annot$Analyte=gsub("@NA$","",paste(newnam,"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
    annot$IsSTD[grep("_ISTD$",newnam)]=TRUE
    annot$IsISO[grep("_ISO$",newnam)]=TRUE
    annot$LevelAnnot[grep("_UNK$",newnam)]=4
    annot$OriginalName=oldnam
    if(!is.null(params$AnnotDB)){
      NewDB=params$AnnotDB
      vnam0=unique(unlist(strsplit(newnam,";")))
      #      print(vnam0)
      lnotfound=unique(vnam0[!vnam0%in%NewDB$GName])
      if(length(lnotfound)>0) cat("Not found in annotation database:\n",lnotfound,"\n",sep=" ")
      l2add=names(NewDB)[names(NewDB)!="GName"]
      toadd=data.frame(sapply(l2add,function(i) sapply(annot$MetName,GRMeta:::.InDBMatchfct,i,NewDB,rmISO=TRUE,rmSTD=FALSE)),stringsAsFactors = F)
      for(i in names(which(sapply(l2add,function(i) is.numeric(NewDB[,i]))))) toadd[,i]=as.numeric(toadd[,i])
      for(i in names(which(sapply(l2add,function(i) is.character(NewDB[,i]))))) toadd[,i]=as.character(toadd[,i])
      annot=cbind(annot,toadd)
      allmat=lapply(allmat,function(x){colnames(x)=annot$Analyte;x})
    }
    rownames(annot)=annot$Analyte
  }
  annot$Method=params$AssayName
  if(params$ordering){
    lso=order(metainfos$sType,fileinfos$Date)
    metainfos=metainfos[lso,]
    fileinfos=fileinfos[lso,]
    allmat=lapply(allmat,function(x) x[lso,,drop=F])
  }
  
  
  l1=c("S/N", "Start Height", "End Height" ,"Start Time", "End Time",  "Del RT")
  l2=c("SNR","Height.Start","Height.End","RT.Start","RT.End","DeltaRT")
  for(i in 1:length(l1)) names(allmat)[names(allmat)==l1[i]]=l2[i]
  
  l2chk=names(allmat)
  if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
  allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})
  
  allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  if(!is.null(ofile)) save(file=ofile,allmat)
  invisible(allmat)
}


