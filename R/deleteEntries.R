deleteEntries<-function(obj,...){
  
  if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
  
  dots <- list(...)
#  print(dots)
  errMsg=paste("Use any of the following for \n","samples:",
               paste(c(names(obj$Meta),names(obj$File)),collapse=" "),"\n","analytes:",
               paste(names(obj$Annot),collapse=" "),"\n",sep="")
  
  
  if(length(dots)==0) stop(errMsg)
 # print(dots)
  for(argname in names(dots)){
    
    if(!argname%in%c(names(obj$Meta),names(obj$File),names(obj$Annot))){
      cat(argname,errMsg,sep=" ")
      next
    }
    
    what=dots[[argname]]
    if(argname%in%c(names(obj$File),names(obj$Meta))){
      if(argname%in%names(obj$Meta)) l2keep=which(!obj$Meta[,argname]%in%what)
      if(argname%in%names(obj$File)) l2keep=which(!obj$File[,argname]%in%what)
      if(length(l2keep)>0 & length(l2keep)<length(obj$Sid)){
        sidori=obj$Sid
        obj$Sid=obj$Sid[l2keep]
        obj$Meta=obj$Meta[l2keep,]
        obj$File=obj$File[l2keep,]
        obj$Data=lapply(obj$Data,function(x) x[l2keep,,drop=F])

        if(!is.null(obj$Eic)) obj$Eic$Sample=obj$Eic$Sample[l2keep,]

        if(any(!sidori%in%obj$Sid)) cat("Samples removed based on ",argname,":\n",sidori[!sidori%in%obj$Sid],"\n",sep=" ")
        if(any(what%in%obj$Sid)) cat("Samples not excluded based on ",argname,":\n",what[what%in%obj$Sid],"\n",sep=" ")
      }
    }
    
    if(argname%in%names(obj$Annot)){
      l2keep=which(!obj$Annot[,argname]%in%what)
      if(length(l2keep)>0 & length(l2keep)<nrow(obj$Annot)){
        anaori=obj$Analyte
        obj$Analyte=obj$Analyte[l2keep]
        obj$Annot=obj$Annot[l2keep,]
        obj$Data=lapply(obj$Data,function(x) x[,l2keep,drop=F])
        
        if(!is.null(obj$Eic)) obj$Eic$File=obj$Eic$File[l2keep,]
        
        if(any(!anaori%in%obj$Analyte)) cat("Analytes removed based on ",argname,":\n",anaori[!anaori%in%obj$Analyte],"\n",sep=" ")
        if(any(what%in%obj$Analyte)) cat("Analytes not excluded based on ",argname,":\n",what[what%in%obj$Analyte],"\n",sep=" ")
      }
    }
  }
  invisible(obj)
}
