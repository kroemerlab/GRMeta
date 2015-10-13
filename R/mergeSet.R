#####################################################################################################
mergeSet<-function(...){
  
  re <- list(...)
  nams=as.character(as.list( match.call() ))[-1]
  tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
  re=re[tokeep]
  names(re)=nams[tokeep]
  if(length(re)==1) return(invisible(re[[1]]))
  
  lmethods=sapply(re,function(x) x$Method)
  #  print(lmethods)
  while(any(table(lmethods)>1)){
    i=names(which(table(lmethods)>1))[1]
    lred=names(which(lmethods==i))
    j=paste(lred,collapse = "")
    cmd=paste(paste(lred,"=re[['",lred,"']]",sep=""),collapse=";")
    eval(parse(text=cmd))
    #      print(lred)
    cmd=paste("re[[lred[1]]]=.mergeDataBatch(",paste(lred,collapse=","),")",sep="")
    eval(parse(text=cmd))
    cmd=paste("rm(list='",paste(lred,collapse="','"),"')",sep="")
    eval(parse(text=cmd))
    re=re[which(!names(re)%in%lred[-1])]
    lmethods=sapply(re,function(x) x$Method)
  }
  #  print(str(re))
  if(length(re)==1) return(invisible(re[[1]]))
  lred=names(re)
  cmd=paste(paste(lred,"=re[['",lred,"']]",sep=""),collapse=";")
  eval(parse(text=cmd))
  cmd=paste("res=.mergeDataMethods(",paste(lred,collapse=","),")",sep="")
  eval(parse(text=cmd))
  cmd=paste("rm(list='",paste(lred,collapse="','"),"')",sep="")
  eval(parse(text=cmd))
  return(invisible(res))  
}


.mergeDataBatch<-function(...){
  
  re <- list(...)
  nams=as.character(as.list( match.call() ))[-1]
  tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
  re=re[tokeep]
  names(re)=nams[tokeep]
  if(length(re)==1) return(re)
  
  cat("************************************\nMerging:", names(re),"\n************************************\n")
  lmethods=unique(sapply(re,function(x) x$Method))
  if(length(lmethods)>1) stop("Applicable to datasets from the same method\n")
  
  cat("Check samples\n")
  lusamp=unlist(lapply(re,function(x) x$Sid))
  lsamp2change=names(which(table(lusamp)>1))
  for(i in names(re)){
    re[[i]]$File$Batch=i
    if(length(lsamp2change)>0)
      re[[i]]=update(re[[i]],what="Sid",formerid =lsamp2change,newid = paste(lsamp2change,i,sep="."),exact = T,swap=F )
  }
  
  metainfos=do.call("rbind",lapply(re,function(x) x$Meta))
  fileinfos=do.call("rbind",lapply(re,function(x) x$File))
  metainfos$InjOrder=order(order(fileinfos$Batch,metainfos$InjOrder))
 
#   lso=order(metainfos$sType,metainfos$InjOrder)
#   metainfos=metainfos[lso,]
#   fileinfos=fileinfos[lso,]
  lusamp=metainfos$Sid
  
  cat("Check analytes\n")
  lumet=unique(unlist(lapply(re,function(x) x$Annot$PutativeName)))
  matexVar=sapply(re,function(x) match(lumet,x$Annot$PutativeName))
  whichds=apply(matexVar,1,function(x) which(!is.na(x))[1])
  annot=re[[1]]$Annot
  for(i in which(whichds>1)) annot[i,]=re[[whichds[i]]]$Annot[matexVar[i,whichds[i]],]
  addAnnot=data.frame(sapply(1:length(re),function(x) re[[x]]$Analyte[matexVar[,x]]),stringsAsFactors=F)
  names(addAnnot)=paste("Analyte.",names(re),sep="")
  annot=cbind(annot,addAnnot)
  
  cat("Merge analytes\n")
  for(i in 1:length(re)) re[[i]]$Data=lapply(re[[i]]$Data,function(x) x[,matexVar[,i]])
  lumeas=unique(unlist(lapply(re,function(x) names(x$Data))))
  allmat=list()
  for(i in lumeas){
    tmp=lapply(re,function(x) x$Data[[i]])
    lnulls=which(sapply(tmp,is.null))
    for(j in lnulls) tmp[[j]]=matrix(NA,nrow=length(re[[j]]$Sid),ncol=length(lumet),dimnames=list(re[[j]]$Sid,lumet))
    allmat[[i]]=do.call("rbind",tmp)
  }
  annot$RT=round(apply(allmat$RT,2,median,na.rm=T),4)
  annot$Analyte=paste(annot$PutativeName,"@",sprintf("%.2f",annot$RT),"-",annot$Method[1],sep="")
  rownames(annot)=annot$Analyte
  allmat=lapply(allmat,function(x){dimnames(x)=list(metainfos$Sid,annot$Analyte);x})
  
  allmat=list(Method=lmethods,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  invisible(allmat)
}

#####################################################################################################
.mergeDataMethods<-function(...){
  
  re <- list(...)
  nams=as.character(as.list( match.call() ))[-1]
  tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
  re=re[tokeep]
  names(re)=nams[tokeep]
  if(length(re)==1) return(re)
  
  cat("************************************\nMerging:", names(re),"\n************************************\n")
  lmethods=sapply(re,function(x) x$Method)
  if(max(table(lmethods))>1) stop("Applicable to datasets from different methods\n")
  names(re)=lmethods
  
  
  lusamp=unique(unlist(lapply(re,function(x) x$Sid)))
  matexSa=sapply(re,function(x) match(lusamp,x$Sid))
  lnmatch=lapply(colnames(matexSa),function(x) lusamp[which(is.na(matexSa[,x]))])
  names(lnmatch)=colnames(matexSa)
  
  if(any(is.na(matexSa))){cat("Missing samples in :\n");print(lnmatch[!sapply(lnmatch,is.null)])}
  whichds=colnames(matexSa)[apply(matexSa,1,function(x) which(!is.na(x))[1])]
  meta=re[[1]]$Meta[,c("Sid","sType")]
  for(i in which(whichds!=colnames(matexSa)[1])) meta[i,]=re[[whichds[i]]]$Meta[matexSa[i,whichds[i]],c("Sid","sType")]
  
  fileinfos=data.frame(Sid=lusamp,stringsAsFactors = FALSE)
  rownames(fileinfos)=lusamp
  for(i in names(re)){
    addfi=re[[i]]$File[,names(re[[i]]$File)!="Sid"]
    addfi$InjOrder=re[[i]]$Meta$InjOrder
    names(addfi)=paste(names(addfi),i,sep=".")
    fileinfos=cbind(fileinfos,addfi[matexSa[,i],])
  }
  neworder=fileinfos[,grep("InjOrder",names(fileinfos))]
  neworder2=apply(neworder,2,function(x) order(order(x)))
  neworder2[is.na(neworder)]=NA
  meta$InjOrder=order(order(rowMeans(neworder2,na.rm=T)))
  
  ###
  lso=order(meta$sType,meta$InjOrder)
  meta=meta[lso,]
  lusamp=lusamp[lso]
  fileinfos=fileinfos[lso,]
  matexSa=sapply(re,function(x) match(lusamp,x$Sid))
  
  ##############
  cat("Check analytes:\n")
  lumet=unlist(lapply(re,function(x) x$Analyte))
  lann=lapply(re,function(x) x$Annot[,grep("^[A-Za-z]+$",names(x$Annot))])
  l2use=unique(unlist(lapply(lann,names)))
  l2usen=names(which(!sapply(lann,function(x) all(l2use%in%names(x)))))
  for(i in l2usen){
    l2add=l2use[!l2use%in%names(lann[[i]])]
    cat(" * adding",l2add,"to annotation of",i,"\n",sep=" ")
    lann[[i]]=cbind(lann[[i]],data.frame(matrix(NA,nrow=nrow(lann[[i]]),ncol=length(l2add),dimnames=list(rownames(lann[[i]]),l2add))))
  }
  annot=do.call("rbind",lann)
  rownames(annot)=lumet=annot$Analyte
  matexVar=sapply(re,function(x) match(lumet,x$Analyte))
  
  ##############
  cat("Merge data:\n")
  #for(i in names(re)) re[[i]]$Data=lapply(re[[i]]$Data,function(x) x[,matexVar[,i]])
  lumeas=unique(unlist(lapply(re,function(x) names(x$Data))))
  allmat=list()
  for(i in lumeas){
    cat(i," ",sep="")
    allmat[[i]]=matrix(NA,nrow=length(lusamp),ncol=length(lumet),dimnames=list(lusamp,lumet))
    for(j in names(re)){
      if(is.null(re[[j]]$Data[[i]]))  cat("\n * not found in ",j,"\n",sep="")
      if(!is.null(re[[j]]$Data[[i]])){
        # cat("    add ",i," from ",j,"\n",sep="")
        lrow=which(!is.na(matexSa[,j]))
        lcol=which(!is.na(matexVar[,j]))
        allmat[[i]][lrow,lcol]=re[[j]]$Data[[i]][matexSa[lrow,j],matexVar[lcol,j]]
      }
    }
  }
  
  allmat=lapply(allmat,function(x){dimnames(x)=list(meta$Sid,annot$Analyte);x})
  
  allmat=list(Method=names(re),Sid=meta$Sid,Analyte=annot$Analyte,Annot=annot,Meta=meta,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  invisible(allmat)
}
