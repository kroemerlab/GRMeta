
loadAgilentDataFlux<-function(ifile,ofile=NULL,params=list()){

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
  ####################################################################################################
  stmets=unique(sort(c(grep("Results$",names(tab)),grep("Results\\.[0-9]+$",names(tab)))))[1]
  mat=apply(as.matrix(tab[-1,stmets:ncol(tab)]),2,function(x) as.numeric(gsub(",",".",x)))
  metnams=colnames(mat)
  lmetinfos=as.matrix(tab[1,])[,stmets:ncol(tab)]
  if(any(grepl(" ",lmetinfos))){
    cat(" ! Spaces in variable names removed !\n")
    lmetinfos=gsub(" ",".",lmetinfos)
  }
  ## 
  lmets=unique(sort(c(grep("Results$",metnams),grep("Results\\.[0-9]+$",metnams))))
  for(i in (1:length(metnams))[-lmets])  metnams[i]=metnams[i-1]
  metnams=gsub("[\\.]{2,}","_",gsub("[\\.]{1,}Results\\.[1-9]$","",gsub("[\\.]{1,}Results$","",metnams)))
  
  niso=q1=q3=rep(NA,length(lmetinfos))
  ismrm=which(grepl("^Qualifier_.*_.*$",metnams))
  isnotmrm=which(grepl("^Qualifier_.*$",metnams))
  isnotmrm=isnotmrm[!isnotmrm%in%ismrm]
  if(length(ismrm)>0){
    q1[ismrm]=as.numeric(sapply(strsplit(metnams[ismrm],"_"),function(x) x[2]))
    q3[ismrm]=as.numeric(sapply(strsplit(metnams[ismrm],"_"),function(x) x[3]))
  }
  if(length(isnotmrm)>0)  q1[isnotmrm]=as.numeric(sapply(strsplit(metnams[isnotmrm],"_"),function(x) x[2]))
  ananams=rep("",length(lmetinfos))
  ananams[!grepl("^Qualifier_.*$",metnams)]=metnams[!grepl("^Qualifier_.*$",metnams)]
  for(i in which(ananams=="")) ananams[i]=ananams[i-1]
  
  for(i in unique(ananams[is.na(q1)])){
    if(all(is.na(q1[which(ananams==i)]))) next
    q1[which(ananams==i & is.na(q1))]=min(q1[which(ananams==i)],na.rm=T)-1
    if(all(is.na(q3[which(ananams==i)]))) next
    q3i=unique(na.omit(q3[which(ananams==i)]))
    q3i=ifelse(length(q3i)==1,q3i,min(q3i)-1)
    q3[which(ananams==i & is.na(q3))]=q3i
  }
  
  for(i in unique(ananams)){
    l=which(ananams==i)
    if(all(is.na(q3[l])) & all(is.na(q1[l]))){niso[l]=0;next}
    if(all(is.na(q3[l]))){niso[l]=round(q1[l]-min(q1[l]));next}
    nq1=round(q1[l]-min(q1[l]))
    nq3=round(q3[l]-min(q3[l]))
    niso[l]=nq1+nq3
  }
  
  metnams2=cleanMetaboNames(metnam=ananams,RegExpr=NA,Syno=NA)$newnam
  newnam=paste(metnams2,"_M",niso,sep="")
  colnames(mat)=paste(newnam,";-;",lmetinfos,sep="")
  
  
  newnamf=factor(newnam,levels=unique(newnam))
  lumetnams=tapply(newnam,newnamf,unique)
  annot=data.frame(Analyte=tapply(newnam,newnamf,unique),MetName=tapply(metnams2,newnamf,unique),
                   IsSTD=FALSE,RT=NA,LevelAnnot=1,
                   OriginalName=tapply(metnams,newnamf,unique),
                   Q1=tapply(q1,newnamf,unique),
                   Q3=tapply(q3,newnamf,unique),
                   Iso=tapply(niso,newnamf,unique),
                   stringsAsFactors = F)
  lso=order(factor(annot$MetName,levels=unique(annot$MetName)),annot$Iso)
  annot=annot[lso,]
  
  #####
  
  ldatatype=sort(unique(lmetinfos))
  allmat=lapply(ldatatype,function(x){
    imat=mat[,which(lmetinfos==x)]
    colnames(imat)=gsub(paste(";-;",x,"$",sep=""),"",colnames(imat))
    imat=imat[,match(annot$Analyte,colnames(imat))]
    dimnames(imat)=list(rownames(metainfos),annot$Analyte)
    imat
  })
  names(allmat)=ldatatype
  
  #########
  rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
  annot$RT=rtmed
  newn=sprintf("%s@%.3f-%s",annot$Analyte,rtmed,params$AssayName)
  annot$Analyte=newn
  NewDB=params$AnnotDB
  vnam0=unique(unlist(strsplit(annot$MetName,";")))
  lnotfound=unique(vnam0[!vnam0%in%NewDB$GName])
  if(length(lnotfound)>0) cat("Not found in annotation database:\n",lnotfound,"\n",sep=" ")
  l2add=names(NewDB)[names(NewDB)!="GName"]
      toadd=data.frame(sapply(l2add,function(i) sapply(annot$MetName,.InDBMatchfct,i,NewDB)),stringsAsFactors = F)
      for(i in names(which(sapply(l2add,function(i) is.numeric(NewDB[,i]))))) toadd[,i]=as.numeric(toadd[,i])
      for(i in names(which(sapply(l2add,function(i) is.character(NewDB[,i]))))) toadd[,i]=as.character(toadd[,i])
      annot=cbind(annot,toadd)
      allmat=lapply(allmat,function(x){colnames(x)=annot$Analyte;x})
    rownames(annot)=annot$Analyte
  annot$Method=params$AssayName
  
 
  if(params$ordering){
    lso=order(metainfos$sType,fileinfos$Date)
    metainfos=metainfos[lso,]
    fileinfos=fileinfos[lso,]
    allmat=lapply(allmat,function(x) x[lso,,drop=F])
  }
  
  l1=c("S/N","BL.End","BL.Start","Int..End","Int..Start")
  l2=c("SNR","Height.Start","Height.End","RT.Start","RT.End")
  for(i in which(l1%in%names(allmat))) names(allmat)[names(allmat)==l1[i]]=l2[i]
  
  l2chk=names(allmat)
  if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
  allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})
  
  m=allmat$RT
  newrt=tapply(annot$RT,annot$MetName,function(x) x[1])
  oldrt=annot$RT
  newrt=newrt[match(annot$MetName,names(newrt))]
  allmat$RT.adjM0=sweep(m,2,oldrt-newrt)
  
  
  allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),c("metaboSet","fluxoSet"))
  if(!is.null(ofile)) save(file=ofile,allmat)
  invisible(allmat)
}



