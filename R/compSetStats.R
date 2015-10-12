
compSetStats<-function(obj,lStats=list(),sType=NULL,Sids=NULL,addtoname=NULL,addtoobject=FALSE){
  
  cvf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*sd(x,na.rm=T)/mean(x,na.rm=T),NA)
  cvrf<-function(x,n=0) ifelse(sum(!is.na(x))>=n,100*mad(x,na.rm=T)/median(x,na.rm=T),NA)
  

  lSa=NULL
  if(!is.null(Sids)) lSa=which(obj$Meta$Sid%in%Sids)
  if(!is.null(sType)) lSa=which(obj$Meta$sType%in%sType)
  if(length(lSa)>0){
    cat(length(lSa),"samples found\n")
    Data=lapply(obj$Data,function(x) x[lSa,,drop=FALSE])
  } else stop("No samples found!\n")
  
  istats=lStats[[1]]
  allst=list()
  for(istats in lStats){
    if(!istats[1]%in%c("Sum","Prop","PropNNA","CV","CVr")){cat("Computation type must be Sum/Prop/PropNNA/CV/CVr for:",istats,"\n");next}
    if(istats[1]=="Sum") eval(parse(text=paste("re<-try(colSums(Data$",istats[2],",na.rm=T),TRUE)",sep="")))
    if(istats[1]=="Prop")  eval(parse(text=paste("re<-try(colSums(Data$",istats[2],",na.rm=T)*100/nrow(!is.na(Data$",istats[2],")),TRUE)",sep="")))
    if(istats[1]=="PropNNA")  eval(parse(text=paste("re<-try(colSums(Data$",istats[2],",na.rm=T)*100/colSums(!is.na(Data$",istats[2],")),TRUE)",sep="")))
    if(istats[1]=="CV") eval(parse(text=paste("re<-try(apply(Data$",istats[2],",2,cvf,",ifelse(!is.na(istats[3]),as.numeric(istats[3]),1),"),TRUE)",sep="")))
    if(istats[1]=="CVr") eval(parse(text=paste("re<-try(apply(Data$",istats[2],",2,cvrf,",ifelse(!is.na(istats[3]),as.numeric(istats[3]),1),"),TRUE)",sep="")))
    if(istats[1]%in%c("Sum","Prop","PropNNA") & !"try-error"%in%class(re))  re[is.na(re)]=0
    if(!"try-error"%in%class(re)) allst[[paste(c(addtoname,istats),collapse="_")]]=re
    if("try-error"%in%class(re)) cat("Cannot compute:",istats,"\n")
  }
  
  res=do.call("cbind",allst)
  if(addtoobject){
    if(!is.null(res)) obj$Annot=cbind(obj$Annot,res)
    invisible(obj)
  } else invisible(res)
    
}

