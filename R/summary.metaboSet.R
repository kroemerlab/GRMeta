
summary.metaboSet<-function(objData){
  
  if(!inherits(objData, "metaboSet")) stop("This is not a metaboSet object")
  
  methods=objData$Method
  #   cat(length(unique(methods)),"methods:",methods,"\n")
  #   
  files=objData$File
  if(any(names(files)=="Date")) names(files)[names(files)=="Date"]=paste(names(files)[names(files)=="Date"],methods,sep=".")
  if(any(names(files)=="Batch")) names(files)[names(files)=="Batch"]=paste(names(files)[names(files)=="Batch"],methods,sep=".")
  alldts=allnum=list()
  stype=objData$Meta$sType
  for(i in methods){
    dts=files[,paste("Date",i,sep=".")]
    ibatch=paste("Batch",i,sep=".")
    if(!ibatch%in%names(file)) cbatch=rep("B1",nrow(files))
    if(ibatch%in%names(file)) cbatch=files[,ibatch]
    alldts[[i]]=tapply(dts,cbatch,range,na.rm=T)
    allnum[[i]]=tapply(stype,cbatch,table)
  }
  tabmet=table(objData$Annot$IsSTD,objData$Annot$Method)
  invisible(list(Meth=methods,Dts=alldts,Distr=allnum,Tab=tabmet))
}
