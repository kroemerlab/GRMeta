update.metaboSet<-function(obj,what="Sid",formerid=NULL,newid=NULL,exact=TRUE,swap=FALSE,verbose=TRUE){
  
  if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
  
  if(!what%in%c("Analyte","Sid")) print("Only Sid or Analyte suported!")
  if(length(formerid)!=length(newid) | is.null(formerid)) stop("former and new inputs are of different length")
  
  oldnames=newnames=obj[[what]]
  if(!swap){
    if(exact) for(i in 1:length(formerid)) newnames[which(oldnames==formerid[i])]=newid[i]
  if(!exact) for(i in 1:length(formerid)) newnames=gsub(formerid[i],newid[i],newnames)
  }
  if(swap) for(i in 1:length(formerid)){
    i1=which(oldnames==formerid[i])
    i2=which(oldnames==newid[i])
    if(length(i1)==1 & length(i2)==1){
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
  
  if(!is.null(obj[["Eic"]])) if(!is.null(obj[["Eic"]]$Samp)) obj[["Eic"]]$Samp$Sid=rownames(obj[["Eic"]]$Samp)=newnames
  
  if(any(newnames!=oldnames)){
   if(verbose){
     cat("Samples with updated id:\n")
     print(rbind(Former=oldnames[newnames!=oldnames],New=newnames[newnames!=oldnames]))
   } else cat("Updating samples id: ",sum(newnames!=oldnames),"/",length(newnames),"\n",sep="")
    
  }
  }
  
  if(what=="Analyte"){
    obj$Analyte=newnames
    obj$Annot$Analyte=rownames(obj$Annot)=newnames
    obj$Data=lapply(obj$Data,function(x){colnames(x)=newnames;x})
    
    if(!is.null(obj[["Eic"]]))  if(!is.null(obj[["Eic"]]$File)) obj[["Eic"]]$File$Analyte=rownames(obj[["Eic"]]$File)=newnames
    if(!is.null(obj[["EicDef"]])) obj[["EicDef"]]$Analyte=rownames(obj[["EicDef"]])=newnames
    
    if(any(newnames!=oldnames)){
      if(verbose){
        cat("Analytes with updated id:\n")
      print(rbind(Former=oldnames[newnames!=oldnames],New=newnames[newnames!=oldnames]))
      } else cat("Updating analytes id: ",sum(newnames!=oldnames),"/",length(newnames),"\n",sep="")
    }
  }
  
  invisible(obj)
}
