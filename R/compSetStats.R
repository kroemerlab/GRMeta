
compSetStats<-function(obj,lStats=list(),whatSamp="Sa",whatQC="QC",add=FALSE){
  
  cvf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*sd(x,na.rm=T)/mean(x,na.rm=T),NA)
  cvrf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*mad(x,na.rm=T)/median(x,na.rm=T),NA)
  

  SaData=obj$Data
  if(length(whatSamp)>0 & whatSamp%in%obj$Meta$Sid) SaData=lapply(obj$Data,function(x) x[which(obj$Meta$Sid%in%whatSamp),,drop=FALSE])
  if(length(whatSamp)==1 & whatSamp%in%obj$Meta$sType) SaData=lapply(obj$Data,function(x) x[which(obj$Meta$sType==whatSamp),,drop=FALSE])
  
  QCData=obj$Data
  if(length(whatQC)>0 & whatQC%in%obj$Meta$Sid) QCData=lapply(obj$Data,function(x) x[which(obj$Meta$Sid%in%whatQC),,drop=FALSE])
  if(length(whatQC)==1 & whatQC%in%obj$Meta$sType) QCData=lapply(obj$Data,function(x) x[which(obj$Meta$sType==whatQC),,drop=FALSE])
  
  istats=lStats[[1]]
  allst=list()
  for(istats in lStats){
    if(!istats[1]%in%c("Sum","Prop","PropNNA","CV","CVr")){cat("Computation type must be Sum/Prop/PropNNA/CV/CVr for:",istats,"\n");next}
    if(!istats[2]%in%c("QC","Sa")){cat("Computation type must be Sa/QC for:",istats,"\n");next}
    
    if(istats[1]=="Sum") eval(parse(text=paste("re<-try(colSums(",istats[2],"Data$",istats[3],",na.rm=T),TRUE)",sep="")))
    if(istats[1]=="Prop")  eval(parse(text=paste("re<-try(colSums(",istats[2],"Data$",istats[3],",na.rm=T)*100/nrow(!is.na(",istats[2],"Data$",istats[3],")),TRUE)",sep="")))
    if(istats[1]=="PropNNA")  eval(parse(text=paste("re<-try(colSums(",istats[2],"Data$",istats[3],",na.rm=T)*100/colSums(!is.na(",istats[2],"Data$",istats[3],")),TRUE)",sep="")))
    if(istats[1]=="CV") eval(parse(text=paste("re<-try(apply(",istats[2],"Data$",istats[3],",2,cvf,",ifelse(!is.na(istats[4]),as.numeric(istats[4]),1),"),TRUE)",sep="")))
    if(istats[1]=="CVr") eval(parse(text=paste("re<-try(apply(",istats[2],"Data$",istats[3],",2,cvrf,",ifelse(!is.na(istats[4]),as.numeric(istats[4]),1),"),TRUE)",sep="")))
    if(istats[1]%in%c("Sum","Prop","PropNNA") & !"try-error"%in%class(re))  re[is.na(re)]=0
    if(!"try-error"%in%class(re)) allst[[paste(istats,collapse="_")]]=re
    if("try-error"%in%class(re)) cat("Cannot compute:",istats,"\n")
  }
  
  res=do.call("cbind",allst)
  if(add){
    if(!is.null(res)) obj$Annot=cbind(obj$Annot,res)
    invisible(obj)
  }

  if(!add) invisible(res)
    
}

