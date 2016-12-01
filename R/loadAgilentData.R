.InDBMatchfct<-function(x,what,NewDB){
  l=which(NewDB$GName%in%strsplit(x,";")[[1]])
  if(length(l)==0) return(NA)
  if(is.numeric(NewDB[,what])) return(mean(NewDB[l,what],na.rm=T))
  l=unique(unlist(strsplit(NewDB[l,what],";"),use.names=FALSE))
  l=l[which(!(l=="" & is.na(l) & l=="NA"))]
  if(length(l)==0) return(NA)
#  if(what=="KEGG") print(c(x,paste(l,collapse=";")))
  paste(l,collapse=";")
  
}



loadAgilentData<-function(ifile,ofile=NULL,params=list()){
  
  
  paramsvals <- paramsParsing()
  if (!missing(params)) paramsvals[names(params)] <- params
  params=paramsvals
  
  tab=read.table(ifile,sep="\t",stringsAsFactors=F,header=T)
  tab=tab[tab[,1]!="",]
  
  fline=strsplit(scan(ifile,sep="\n",what="raw"),"\t")[[1]]
  names(tab)[which(fline!="")]=fline[which(fline!="")]
  
  #########################################################################################################################
  whichdate=which(as.matrix(tab[1,])[1,]==params$TimeCol )
  
  dts0=(as.matrix(tab[,whichdate]))[-1,1]
  
  dts=sapply(strsplit(dts0," "),function(x) x[1])
  hrs=sapply(strsplit(dts0," "),function(x) paste(x[-1],collapse=" "))
  hrs=substr(strptime(hrs,"%I:%M %p"),12,19)
  if(all(is.na(hrs))) hrs=paste(gsub(".* ","",dts0),":00",sep="")
  dts=chron(dates.=dts,times.=hrs,format=c(dates="m/d/y",times="h:m:s"))
  
  #########################################################################################################################
  
  whichfile=which(as.matrix(tab[1,])[1,]==params$FileCol )
  
  filenam=(as.matrix(tab[,whichfile]))[-1,1]
  
  if(!is.null(params$NameCol)){
    whichname= which(as.matrix(tab[1,])[1,]==params$NameCol )
    if(length(whichname)==0) params$NameCol=NULL else  nams=(as.matrix(tab[,whichname]))[-1,1]
  } 
  if(is.null(params$NameCol))  nams=gsub("\\.[dD]$","",filenam)
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
  fileinfos=data.frame(File=filenam,Date=dts,Name=nams,Sid=sid,stringsAsFactors=FALSE)
  rownames(fileinfos)=fileinfos$Sid
  if(!is.null(params$Batch)) fileinfos$Batch=params$Batch
  
  lmets=grep("Results$",names(tab))
  lmetinfos=as.matrix(tab[1,])[,lmets[1]:ncol(tab)]
  metnams=rep("",length(lmetinfos))
  metnams[lmets-lmets[1]+1]=gsub(" Results$","",names(tab)[lmets])
  for(i in which(metnams=="")) metnams[i]=metnams[i-1]
  lumetnams=unique(metnams)
  
  mat=apply(as.matrix(tab[-1,lmets[1]:ncol(tab)]),2,
            function(x) as.numeric(gsub(",",".",gsub("Infinity","Inf",x))))
  
  ldatatype=sort(unique(lmetinfos))
  allmat=lapply(ldatatype,function(x){
    imat=mat[,which(lmetinfos==x)]
    imat=imat[,match(lumetnams,metnams[which(lmetinfos==x)])]
    dimnames(imat)=list(rownames(metainfos),lumetnams)
    imat
  })
  names(allmat)=ldatatype
  
  #########
  rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
  nm=gsub("@NA$","",paste(lumetnams,"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
  annot=data.frame(Analyte=nm,MetName=lumetnams,IsSTD=FALSE,RT=rtmed,LevelAnnot=1,stringsAsFactors = F)
  newnam=oldnam=lumetnams
  
  if(params$checkNams){
    newnam=cleanMetaboNames(oldnam,RegExpr = NA,Syno = NA)$newnam
    annot$MetName=newnam
    annot$Analyte=gsub("@NA$","",paste(newnam,"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
    annot$IsSTD[grep("_ISTD",newnam)]=TRUE
    annot$OriginalName=oldnam
    if(!is.null(params$AnnotDB)){
      NewDB=params$AnnotDB
      vnam0=unique(unlist(strsplit(newnam,";")))
#      print(vnam0)
      lnotfound=unique(vnam0[!vnam0%in%NewDB$GName])
      if(length(lnotfound)>0) cat("Not found in annotation database:\n",lnotfound,"\n",sep=" ")
      l2add=names(NewDB)[names(NewDB)!="GName"]
      toadd=data.frame(sapply(l2add,function(i) sapply(annot$MetName,.InDBMatchfct,i,NewDB)),stringsAsFactors = F)
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
  
  l1=c("S/N","BL Start","BL End","Int. Start","Int. End")
  l2=c("SNR","Height.Start","Height.End","RT.Start","RT.End")
  for(i in 1:length(l1)) names(allmat)[names(allmat)==l1[i]]=l2[i]
  
  l2chk=names(allmat)
  if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
  allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})
  
  allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  if(!is.null(ofile)) save(file=ofile,allmat)
  invisible(allmat)
}





.loadAgilentData<-function(ifile,ofile=NULL,params=list()){
  

  paramsvals <- paramsParsing()
  if (!missing(params)) paramsvals[names(params)] <- params
  params=paramsvals
  
  tab=read.table(ifile,sep="\t",stringsAsFactors=F,header=T)#[1:74,]
  tab=tab[tab[,1]!="",]
  
  #########################################################################################################################
  whichdate=which(as.matrix(tab[1,])[1,]==params$TimeCol )
  
  dts0=(as.matrix(tab[,whichdate]))[-1,1]
  
  dts=sapply(strsplit(dts0," "),function(x) x[1])
  hrs=sapply(strsplit(dts0," "),function(x) paste(x[-1],collapse=" "))
  hrs=substr(strptime(hrs,"%I:%M %p"),12,19)
  if(all(is.na(hrs))) hrs=paste(gsub(".* ","",dts0),":00",sep="")
  dts=chron(dates.=dts,times.=hrs,format=c(dates="m/d/y",times="h:m:s"))
  
  #########################################################################################################################
  
  whichfile=which(as.matrix(tab[1,])[1,]==params$FileCol )
  
  filenam=(as.matrix(tab[,whichfile]))[-1,1]
  
  if(!is.null(params$NameCol)){
    whichname= which(as.matrix(tab[1,])[1,]==params$NameCol )
    if(length(whichname)==0) params$NameCol=NULL else  nams=(as.matrix(tab[,whichname]))[-1,1]
  } 
  if(is.null(params$NameCol))  nams=gsub("\\.[dD]$","",filenam)
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
  fileinfos=data.frame(File=filenam,Date=dts,Name=nams,Sid=sid,stringsAsFactors=FALSE)
  rownames(fileinfos)=fileinfos$Sid
  if(!is.null(params$Batch)) fileinfos$Batch=params$Batch

  lmets=grep("Results$",names(tab))
  lmetinfos=as.matrix(tab[1,])[,lmets[1]:ncol(tab)]
  metnams=rep("",length(lmetinfos))
  metnams[lmets-lmets[1]+1]=gsub("\\.Results$","",names(tab)[lmets])
  for(i in which(metnams=="")) metnams[i]=metnams[i-1]
  lumetnams=unique(metnams)
  
  mat=apply(as.matrix(tab[-1,lmets[1]:ncol(tab)]),2,function(x) as.numeric(gsub(",",".",x)))
  
  ldatatype=sort(unique(lmetinfos))
  allmat=lapply(ldatatype,function(x){
    imat=mat[,which(lmetinfos==x)]
    imat=imat[,match(lumetnams,metnams[which(lmetinfos==x)])]
    dimnames(imat)=list(rownames(metainfos),lumetnams)
    imat
  })
  names(allmat)=ldatatype
  
  #########
  rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
  nm=gsub("@NA$","",paste(lumetnams,"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
  annot=data.frame(Analyte=nm,MetName=lumetnams,IsSTD=FALSE,RT=rtmed,LevelAnnot=1,stringsAsFactors = F)
  newnam=oldnam=lumetnams
  
  if(params$checkNams){
    newnam=cleanMetaboNames(oldnam,RegExpr = NA,Syno = NA)$newnam
    annot$MetName=newnam
    annot$Analyte=gsub("@NA$","",paste(newnam,"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
    annot$IsSTD[grep("_ISTD",newnam)]=TRUE
    annot$OriginalName=oldnam
    if(!is.null(params$AnnotDB)){
      NewDB=params$AnnotDB
      vnam0=unique(unlist(strsplit(newnam,";")))
      lnotfound=unique(vnam0[!vnam0%in%NewDB$GName])
      if(length(lnotfound)>0) cat("Not found in annotation database:\n",lnotfound,"\n",sep=" ")
      l2add=names(NewDB)[names(NewDB)!="GName"]
      toadd=data.frame(sapply(l2add,function(i) sapply(annot$MetName,.InDBMatchfct,i,NewDB)),stringsAsFactors = F)
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
  
  l1=c("S/N","BL Start","BL End","Int. Start","Int. End")
  l2=c("SNR","Height.Start","Height.End","RT.Start","RT.End")
  for(i in 1:length(l1)) names(allmat)[names(allmat)==l1[i]]=l2[i]
  
  l2chk=names(allmat)
  if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
  allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})
  
  allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  if(!is.null(ofile)) save(file=ofile,allmat)
  invisible(allmat)
}


