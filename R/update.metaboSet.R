update.metaboSet<-function(obj,what="Sid",formerid=NULL,newid=NULL,exact=TRUE,swap=FALSE){
  
  if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
  
  if(!what%in%c("Analyte","Sid")) print("Only Sid or Analyte suported!")
  if(length(formerid)!=length(newid) | is.null(formerid)) stop("former and new inputs are of different length")
  
  oldnames=newnames=obj[[what]]
  if(exact & !swap) for(i in 1:length(formerid)) newnames[which(oldnames==formerid[i])]=newid[i]
  if(!exact & !swap) for(i in 1:length(formerid)) newnames=gsub(formerid[i],newid[i],newnames)
  if(swap) for(i in 1:length(formerid)){
    i1=which(oldnames==formerid[i])
    i2=which(oldnames==newid[i])
    if(length(i1)>0 & length(i2)){
      newnames[i1]=newid[i]
      newnames[i2]=formerid[i]
    }
  }
  if(any(table(newnames)>1)) stop(paste(c("duplicated names: ",names(which(table(newnames)>1))), collapse=" "))
  
  if(what=="Sid"){
  obj$Sid=newnames
  obj$Meta$Sid=rownames(obj$Meta)=newnames
  obj$File$Sid=rownames(obj$File)=newnames
  obj$Data=lapply(obj$Data,function(x){rownames(x)=newnames;x})
  if(any(newnames!=oldnames)) print(rbind(Old=oldnames[newnames!=oldnames],New=newnames[newnames!=oldnames]))
  }
  
  if(what=="Analyte"){
    obj$Analyte=newnames
    obj$Annot$Analyte=rownames(obj$Annot)=newnames
    obj$Data=lapply(obj$Data,function(x){colnames(x)=newnames;x})
    if(any(newnames!=oldnames)) print(rbind(Old=oldnames[newnames!=oldnames],New=newnames[newnames!=oldnames]))
  }
  
  invisible(obj)
}
