sortEntries<-function(obj,...){
  
  if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
  
  dots <- list(...)

  errMsg=paste("Use any of the following for \n","samples:",
               paste(c(names(obj$Meta),names(obj$File)),collapse=" "),"\n","analytes:",
               paste(names(obj$Annot),collapse=" "),"\n",sep="")
  
  
  if(length(dots)==0) stop(errMsg)

  ssort=asort=NULL
  
  for(argname in names(dots)){
    
    if(!argname%in%c(names(obj$Meta),names(obj$File),names(obj$Annot))){
      cat(argname,errMsg,sep=" ")
      next
    }
    
    what=dots[[argname]]
    if(!is.logical(what)) what=TRUE
    
    if(argname%in%c(names(obj$File),names(obj$Meta))){
      if(argname%in%names(obj$Meta)) ssort=cbind(ssort,xtfrm(obj$Meta[,argname])*ifelse(what,1,-1))
      if(argname%in%names(obj$File)) ssort=cbind(ssort,xtfrm(obj$File[,argname])*ifelse(what,1,-1))
      colnames(ssort)[ncol(ssort)]=paste(argname,":",what,sep="")
    }
    
    if(argname%in%names(obj$Annot)){
      asort=cbind(asort,xtfrm(obj$Annot[,argname])*ifelse(what,1,-1))
      colnames(asort)[ncol(asort)]=paste(argname,":",what,sep="")
    }
  }
  
  
  if(!is.null(ssort)){
    lso=do.call(order, as.data.frame(ssort))
    obj$Sid=obj$Sid[lso]
    obj$Meta=obj$Meta[lso,]
    obj$File=obj$File[lso,]
    obj$Data=lapply(obj$Data,function(x) x[lso,,drop=F])
    if(!is.null(obj$Eic)) obj$Eic$Sample=obj$Eic$Sample[lso,]
    cat("Sorting samples based on:",colnames(ssort),"\n",sep=" ")
  }  
  if(!is.null(asort)){
    lso=do.call(order, as.data.frame(asort))
    obj$Analyte=obj$Analyte[lso]
    obj$Annot=obj$Annot[lso,]
    obj$Data=lapply(obj$Data,function(x) x[,lso,drop=F])
    if(!is.null(obj$Eic)) obj$Eic$File=obj$Eic$File[lso,]
    cat("Sorting analytes based on",colnames(asort),"\n",sep=" ")
  }
  
  invisible(obj)
}
