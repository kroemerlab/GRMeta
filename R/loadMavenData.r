.GRgetcomptime<-function(ifile){
  if(!file.exists(ifile)) return(NA)
  infile<-try(suppressWarnings(suppressMessages(readLines(ifile))),TRUE)
  if("try-error"%in%class(infile)) return(NA)
  dts=grep('completionTime',infile)
  if(length(dts)==0) return(NA)
  dts=infile[dts[1]]
  rm(list="infile")
  dts=strsplit(gsub('\"','',regmatches(dts,regexpr('\"(.*)\"',dts))),"T")[[1]]
  chron(dts[1],dts[2],format=c(dates="Y-M-D",times="h:m:s"))
}

loadMavenData<-function(ifile,ofile=NULL,stdData=NULL,datadir=NA,params=list()){
  
  ####
  infct<-function(cefi){
    require(XML)
    system.time(doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))[[2]])
    dproc=xmlChildren(xmlChildren(doc)$dataProcessing)
    dts=xmlAttrs(dproc$software)["completionTime"]
    dts = strsplit(gsub("\"", "", dts), "T")[[1]]
    unname(chron(dts[1], dts[2], format = c(dates = "Y-M-D", times = "h:m:s")))
  }
  ####
  
  paramsvals <-paramsParsing()
  if (!missing(params)) paramsvals[names(params)] <- params
  params=paramsvals
  
  
  dat <- xmlRoot(xmlTreeParse(ifile))
  
  #####################
  # samples

  sampids=t(xmlSApply(dat[["samples"]],xmlAttrs))
  
  filenam=sampids[,"filename"]
  sampnams=sampids[,"name"]
  
  dts=rep(NA,length(sampnams))
  ### check if RAW files from Thermo
  if(is.character(datadir)) if(dir.exists(datadir) & any(is.na(dts))) for(ftyp in c(".RAW",".raw",".Raw")){
  lfiles=file.path(datadir, paste0(sampnams,ftyp))
  if(!any(is.na(dts) & file.exists(lfiles)) ) next
    l=which(is.na(dts) & file.exists(lfiles))
    dts[l]=as.chron(file.info(lfiles[l])[,"mtime"], format = c(dates = "Y-M-D", times = "h:m:s"))
  }
  ### check if mzdata from Agilent
  if(is.character(datadir)) if(dir.exists(datadir) & any(is.na(dts))){
    lfiles=file.path(datadir, paste0(sampnams,".mzdata.xml"))
    if(any(is.na(dts) & file.exists(lfiles)) ){
    l=which(is.na(dts) & file.exists(lfiles))
    dts[l]=sapply(lfiles[l],infct)
    }
  }
  if(any(!is.na(dts))) dts=as.chron(dts)
  
  ########################################
  nams=gsub("\\.[dD]$","",sampnams)
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
  
  metainfos=data.frame(Sid=sid,sType=styp,InjOrder=as.numeric(sampids[,"sampleOrder"]),stringsAsFactors=FALSE)
  rownames(metainfos)=metainfos$Sid
  fileinfos=data.frame(File=filenam,Date=dts,Name=nams,Sid=sid,stringsAsFactors=FALSE)
  rownames(fileinfos)=fileinfos$Sid
  if(!is.null(params$Batch)) fileinfos$Batch=params$Batch
  
  #####################
  pkids= t(xmlSApply(dat[["PeakGroups"]],xmlAttrs))
  
  #################
  ## Make dat
  m0=matrix(NA,nrow=nrow(sampids),ncol=nrow(pkids),dimnames=list(unname(sampnams),unname(pkids[,1])))
  lnames=c( "rt","rtmin", "rtmax","medianMz","mzmin","mzmax","peakArea" , "peakIntensity")
  lnames2=c( "RT","RT.min", "RT.max","MZ","MZ.min","MZ.max","Area" , "Height")
  allmat=lapply(lnames,function(x) m0)
  names(allmat)=lnames
  lgrps=xmlChildren(dat[["PeakGroups"]])
  
  cat("Found:",length(lgrps),"peak groups ")
  for(k in 1:length(lgrps)){
    if(k%%100==0) cat("x")
    igrp=xmlAttrs(lgrps[[k]])["groupId"]
    mgrp=xmlSApply(lgrps[[k]],xmlAttrs)
    for(i in rownames(mgrp)[rownames(mgrp)%in%lnames]) allmat[[i]][mgrp["sample",],igrp]=as.numeric(mgrp[i,])
  }
  names(allmat)=lnames2
  cat("\n")
  
  #################
  ## Make annot
  ## double check redundant analyte names
  .rendups<-function(oldnames){
    l=names(which(table(oldnames)>1))
    if(length(l)==0) return(oldnames)
    cat("  ! duplicated peaknames:")
    for(i in l){
      li=which(oldnames==i)
      cat( " ",oldnames[li][1],sep="")
      oldnames[li]=paste(gsub("(.*)@(.*)","\\1",oldnames[li]),"_D",0:(length(li)-1),
                         gsub("(.*)@(.*)","@\\2",oldnames[li]),sep="")
    }
    cat("\n")
    return(oldnames)
  } 
  
  
  rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
  mzmed=round(apply(allmat$MZ,2,median,na.rm=T),5)
  if(!"compoundName"%in%colnames(pkids)){
    nm=sprintf("%.4f@%.3f-%s",mzmed,rtmed,params$AssayName)
    nm=.rendups(nm)
    annot=data.frame(Analyte=nm,MetName=NA,IsSTD=FALSE,IsISO=FALSE,RT=rtmed,MZ=mzmed,LevelAnnot=4,Method=params$AssayName,stringsAsFactors = F)
  }
  
  ##############
  # do annot here
  if("compoundName"%in%colnames(pkids)) {
  nm=gsub("@NA$","",paste(unname(pkids[,"compoundName"]),"@",sprintf("%.3f",rtmed),"-",params$AssayName,sep=""))
  nm=.rendups(nm)
  if(!is.null(stdData)) annot=data.frame(Analyte=nm,MetName=unname(pkids[,"compoundName"]),IsSTD=FALSE,IsISO=FALSE,
                                    RT=rtmed,DRT=NA,MZ=mzmed,DPPM=NA,LevelAnnot=1,Method=params$AssayName,stringsAsFactors = F)
  if(is.null(stdData)) annot=data.frame(Analyte=nm,MetName=unname(pkids[,"compoundName"]),IsSTD=FALSE,IsISO=FALSE,
                                    RT=rtmed,MZ=mzmed,LevelAnnot=1,Method=params$AssayName,stringsAsFactors = F)
  newnam=oldnam=unname(pkids[,"compoundName"])
  
  if(params$checkNams){
    newnam=cleanMetaboNames(oldnam,RegExpr = NA,Syno = NA)$newnam
    annot$MetName=newnam
    nm=gsub("@NA$","",paste(newnam,"@",sprintf("%.3f",rtmed),"-",params$AssayName,sep=""))
    nm=.rendups(nm)
    annot$Analyte=nm
    annot$IsSTD[grep("_ISTD",newnam)]=TRUE
    annot$IsISO[grep("_ISO",newnam)]=TRUE
    annot$OriginalName=oldnam
    if(!is.null(params$AnnotDB)){
      NewDB=params$AnnotDB
      vnam0=unique(unlist(strsplit(newnam,";")))
      lnotfound=unique(vnam0[!vnam0%in%NewDB$GName])
      if(length(lnotfound)>0) cat("Not found in annotation database:\n",lnotfound,"\n",sep=" ")
      l2add=names(NewDB)[names(NewDB)!="GName"]
      toadd=data.frame(sapply(l2add,function(i) sapply(annot$MetName,GRMeta:::.InDBMatchfct,i,NewDB,rmISO=TRUE,rmSTD=FALSE)),stringsAsFactors = F)
      for(i in names(which(sapply(l2add,function(i) is.numeric(NewDB[,i]))))) toadd[,i]=as.numeric(toadd[,i])
      for(i in names(which(sapply(l2add,function(i) is.character(NewDB[,i]))))) toadd[,i]=as.character(toadd[,i])
      annot=cbind(annot,toadd)
    }
  }
  
  if(!is.null(stdData)){
    for(imet in annot$MetName[which(annot$MetName%in%stdData$GName)]){
      tmp=stdData[which(stdData$GName==imet),]
      l=which(annot$MetName==imet)
      annot$DRT[l]=round(annot$RT[l]-median(tmp$RT),3)
      annot$DPPM[l]=round(10^6*(annot$MZ[l]-median(tmp$IP))/annot$MZ[l],2)
    }
  }
  
  }
  
  allmat=lapply(allmat,function(x){colnames(x)=annot$Analyte;x})
  rownames(annot)=annot$Analyte
  
  
  if(params$ordering){
    lso=order(metainfos$sType,fileinfos$Date)
    metainfos=metainfos[lso,]
    fileinfos=fileinfos[lso,]
    allmat=lapply(allmat,function(x) x[lso,,drop=F])
  }
  
  l2chk=names(allmat)
  if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
  allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})
  
  allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  


if(!is.null(ofile)) save(file=ofile,allmat)
invisible(allmat)
}



