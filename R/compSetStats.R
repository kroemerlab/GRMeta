
compSetStats<-function(obj,lStats=list(),whatSamp="Sa",whatQC="QC",add=FALSE){
  
  cvf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*sd(x,na.rm=T)/mean(x,na.rm=T),NA)
  cvrf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*mad(x,na.rm=T)/median(x,na.rm=T),NA)
  

  lSampType=lSa=lQC=NULL
  if(length(whatSamp)>0 & any(whatSamp%in%obj$Meta$Sid)) lSa=which(obj$Meta$Sid%in%whatSamp)
  if(length(whatSamp)==1 & whatSamp[1]%in%obj$Meta$sType)  lSa=which(obj$Meta$sType==whatSamp)
  if(length(lSa)>0){
    lSampType=c(lSampType,"Sa")
    cat(length(lSa),"sample found\n")
    SaData=lapply(obj$Data,function(x) x[lSa,,drop=FALSE])
  }
  
  
  if(length(whatQC)>0 & any(whatQC%in%obj$Meta$Sid)) lQC=which(obj$Meta$Sid%in%whatQC)
  if(length(whatQC)==1 & whatQC[1]%in%obj$Meta$sType)  lQC=which(obj$Meta$sType==whatQC)
  if(length(lQC)>0){
    lSampType=c(lSampType,"QC")
    cat(length(lQC),"QC found\n")
    QCData=lapply(obj$Data,function(x) x[lQC,,drop=FALSE])
  }
  if(length(lSampType)==0) stop("No samples found!\n")
  istats=lStats[[1]]
  allst=list()
  for(istats in lStats){
    if(!istats[1]%in%c("Sum","Prop","PropNNA","CV","CVr")){cat("Computation type must be Sum/Prop/PropNNA/CV/CVr for:",istats,"\n");next}
    if(!istats[2]%in%lSampType){cat("Computation type must be ",paste(lSampType,collapse="/")," for:",istats,"\n");next}
    
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
  } else invisible(res)
    
}

