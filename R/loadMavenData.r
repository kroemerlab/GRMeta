# 
# library(XML)
# library(GRMeta)
# 
# ifile="./exampleQTOF.mzPeaks.mzroll"

loadMavenData<-function(ifile,ofile=NULL,params=list()){
  
  require(XML)
params=list()
paramsvals <- GRMeta:::paramsParsing()
if (!missing(params)) paramsvals[names(params)] <- params
params=paramsvals


data <- xmlRoot(xmlTreeParse(ifile))

#####################
# samples
sampids= t(xmlSApply(data[["samples"]],xmlAttrs))

filenam=sampids[,"filename"]
nams=gsub("\\.[dD]$","",sampids[,"name"])
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
fileinfos=data.frame(File=filenam,Date=NA,Name=nams,Sid=sid,stringsAsFactors=FALSE)
rownames(fileinfos)=fileinfos$Sid
if(!is.null(params$Batch)) fileinfos$Batch=params$Batch

#####################
pkids= t(xmlSApply(data[["PeakGroups"]],xmlAttrs))

#################
## Make data
m0=matrix(NA,nrow=length(snames),ncol=nrow(pkids),dimnames=list(unname(sampids[,"name"]),unname(pkids[,1])))
lnames=c( "rt","rtmin", "rtmax","mzmin","mzmax","peakArea" , "peakIntensity")
lnames2=c( "RT","RT.min", "RT.max","MZ.min","MZ.max","Area" , "Height")
allmat=lapply(lnames,function(x) m0)
names(allmat)=lnames
lgrps=xmlChildren(data[["PeakGroups"]])

for(k in 1:length(lgrps)){
  cat(".")
  igrp=xmlAttrs(lgrps[[k]])["groupId"]
mgrp=sapply(lgrps[[k]],xmlAttrs)
for(i in rownames(mgrp)[rownames(mgrp)%in%lnames]) allmat[[i]][mgrp["sample",],igrp]=as.numeric(mgrp[i,])
}
names(allmat)=lnames2

#################
## Make annot

rtmed=round(apply(allmat$RT,2,median,na.rm=T),4)
nm=gsub("@NA$","",paste(unname(pkids[,"compoundName"]),"@",sprintf("%.2f",rtmed),"-",params$AssayName,sep=""))
annot=data.frame(Analyte=nm,MetName=unname(pkids[,"compoundName"]),IsSTD=FALSE,RT=rtmed,LevelAnnot=1,stringsAsFactors = F)
newnam=oldnam=unname(pkids[,"compoundName"])

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

l2chk=names(allmat)
if(is.null(params$nozeroscheck)) l2chk=l2chk[!l2chk%in%params$nozeroscheck]
allmat[l2chk]=lapply(allmat[l2chk],function(x){x[which(x<=0)]=NA;x})

allmat=list(Method=params$AssayName,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
class(allmat)=append(class(allmat),"metaboSet")
if(!is.null(ofile)) save(file=ofile,allmat)
invisible(allmat)
}



