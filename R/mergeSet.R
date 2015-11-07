
### change the .merge to importe list t avoid call mergetSet(re,re1) -> 


#####################################################################################################
mergeSet<-function(...){
  
  re <- list(...)
  nams=as.character(as.list( match.call() ))[-1]
  tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
  re=re[tokeep]
  names(re)=nams[tokeep]
  if(length(re)==1) return(invisible(re[[1]]))
  
  lmethods=sapply(re,function(x) x$Method)
  while(any(table(lmethods)>1)){
    i=names(which(table(lmethods)>1))[1]
    lred=names(which(lmethods==i))
    j=paste(lred,collapse = "")
          print(lred)
    re[[lred[1]]]=.mergeDataBatch(re[lred])
    re=re[which(!names(re)%in%lred[-1])]
    lmethods=sapply(re,function(x) x$Method)
  }
  #  print(str(re))
  if(length(re)==1) return(invisible(re[[1]]))
  lred=names(re)
  res=.mergeDataMethods(re[lred])
  return(invisible(res))  
}
######################################################
.joinDF<-function(ldf,lmin){
  ldfn=unique(unlist(lapply(ldf,names)))
  ldfn=c(lmin,ldfn[!ldfn%in%lmin])
  for(i in 1:length(ldf)){
    lmiss=ldfn[!ldfn%in%names(ldf[[i]])]
    if(length(lmiss)>0) 
      ldf[[i]]=cbind(ldf[[i]],matrix(NA,nrow=nrow(ldf[[i]]),ncol=length(lmiss),dimnames=list(rownames(ldf[[i]]),lmiss)))
    ldf[[i]]=ldf[[i]][,ldfn]
  }
  ndf=do.call("rbind",ldf)
  chkfac=sapply(ldf,function(y) sapply(names(y),function(x) is.factor(y[,x])))
  ndf
}

.joinDF2<-function(ldf,lmin,matex){
  ldfn=unique(unlist(lapply(ldf,names)))
  ldfn=c(lmin,ldfn[!ldfn%in%lmin])
  for(i in 1:length(ldf)){
    lmiss=ldfn[!ldfn%in%names(ldf[[i]])]
    if(length(lmiss)>0) 
      ldf[[i]]=cbind(ldf[[i]],matrix(NA,nrow=nrow(ldf[[i]]),ncol=length(lmiss),dimnames=list(rownames(ldf[[i]]),lmiss)))
    
    ldf[[i]]=ldf[[i]][matex[,i],ldfn]
  }
  fdf=ldf[[1]]
  for(j in 2:length(ldf)){
    torep=which(is.na(fdf) & ! is.na(ldf[[j]]),arr.ind =TRUE)
    if(nrow(torep)>0)  for(k in 1:nrow(torep))
      fdf[torep[k,1],torep[k,2]]=ldf[[j]][torep[k,1],torep[k,2]]
  }
  fdf
}
#################################################################

.mergeDataBatch<-function(re){
  
  ## take a list a input -> better than names
#   
#   re <- list(...)
#   nams=as.character(as.list( match.call() ))[-1]
#   tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
#   re=re[tokeep]
#   names(re)=nams[tokeep]
  if(length(re)==1) return(re)
  
  cat("************************************\nMerging:", names(re),"\n************************************\n")
  lmethods=unique(sapply(re,function(x) x$Method))
  if(length(lmethods)>1) stop("Applicable to datasets from the same method\n")
  
  cat("Creating batches:")
  lbatch=c()
  for(i in names(re)){
    if(!is.null(re[[i]]$File$Batch)){
      if(any(re[[i]]$File$Batch%in%lbatch)){
        l=re[[i]]$File$Batch%in%lbatch
        re[[i]]$File$Batch[l]=paste(re[[i]]$File$Batch[l],".1",sep="")
      }
    }
    if(is.null(re[[i]]$File$Batch))
      re[[i]]$File$Batch=paste(i,ifelse(i%in%lbatch,".1",""),sep="")

    lbatch=append(lbatch,unique(re[[i]]$File$Batch))
  }
  cat(lbatch,"\n")
  
  cat("Updating sample ids\n")
  lusamp=unlist(lapply(re,function(x) x$Sid))
  lsamp2change=names(which(table(lusamp)>1))
  for(i in names(re)){
    if(any(re[[i]]$Sid%in%lsamp2change)){
      l=which(re[[i]]$Sid%in%lsamp2change)
      re[[i]]=update(re[[i]],what="Sid",formerid =re[[i]]$Sid[l],
                newid = paste(re[[i]]$Sid[l],re[[i]]$File$Batch[l],sep="."),exact = T,swap=F )
    }
  }
  
  cat("Merging sample ids\n")
  metainfos=.joinDF(lapply(re,function(x) x$Meta),c("Sid","sType","InjOrder"))
  fileinfos=.joinDF(lapply(re,function(x) x$File),c( "File" ,"Date",  "Name",  "Sid", "Batch"))
  metainfos$InjOrder=order(order(fileinfos$Date,factor(fileinfos$Batch,levels=lbatch),metainfos$InjOrder))
  rownames(fileinfos)=rownames(metainfos)=metainfos$Sid
  lusamp=metainfos$Sid
  
#   lso=order(metainfos$sType,metainfos$InjOrder)
#   metainfos=metainfos[lso,]
#   fileinfos=fileinfos[lso,]
  
  lvar2join= c("Analyte","MetName","IsSTD","RT", "LevelAnnot","OriginalName")
  if("fluxoSet" %in%class(re[[1]])) lvar2join= c("Analyte","MetName","IsSTD","RT", "Iso","LevelAnnot","OriginalName")
  
  if(!"fluxoSet" %in%class(re[[1]])) lvarfct=function(x) x$Annot$MetName
  if("fluxoSet" %in%class(re[[1]]))  lvarfct=function(x) paste(x$Annot$MetName,x$Annot$Iso,sep="_M")
  lumet=unique(unlist(lapply(re,lvarfct)))
  cat("Found",length(lumet),"unique analytes\n")
  matexVar=sapply(re,function(x) match(lumet,lvarfct(x)))
  for(i in colnames(matexVar)) if(any(is.na(matexVar[,i]))) cat("  *missing in",i,":",lumet[which(is.na(matexVar[,i]))],"\n")
  annot=.joinDF2(lapply(re,function(x) x$Annot),lvar2join,matexVar)
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
  if(!"fluxoSet" %in%class(re[[1]]))  annot$Analyte=paste(annot$MetName,"@",sprintf("%.2f",annot$RT),"-",annot$Method[1],sep="")
  if("fluxoSet" %in%class(re[[1]]))  annot$Analyte=paste(annot$MetName,"_M",annot$Iso,"@",sprintf("%.2f",annot$RT),"-",annot$Method[1],sep="")
  rownames(annot)=annot$Analyte
  allmat=lapply(allmat,function(x){dimnames(x)=list(metainfos$Sid,annot$Analyte);x})
  
  allmat=list(Method=lmethods,Sid=metainfos$Sid,Analyte=annot$Analyte,Annot=annot,Meta=metainfos,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  if("fluxoSet"%in%class(re[[1]])) class(allmat)=append(class(allmat),"fluxoSet")
  invisible(allmat)
}


#####################################################################################################
.mergeDataMethods<-function(re){
  
  ## take a list a input -> better than names
#  re <- list(...)
 # nams=as.character(as.list( match.call() ))[-1]
#  tokeep=sapply(unique(nams),function(x) which(nams==x)[1])
#  re=re[tokeep]
 # names(re)=nams[tokeep]
  if(length(re)==1) return(re)
  
  cat("************************************\nMerging:", names(re),"\n************************************\n")
  lmethods=sapply(re,function(x) x$Method)
  if(max(table(lmethods))>1) stop("Applicable to datasets from different methods\n")
  names(re)=lmethods
  
  cat("Checking samples\n")
  lusamp=unique(unlist(lapply(re,function(x) x$Sid)))
  matexSa=sapply(re,function(x) match(lusamp,x$Sid))
  ##
  lnmatch=lapply(colnames(matexSa),function(x) lusamp[which(is.na(matexSa[,x]))])
  names(lnmatch)=colnames(matexSa)
  if(any(is.na(matexSa))){cat("Missing samples in :\n");print(lnmatch[!sapply(lnmatch,is.null)])}
  if(all(rowSums(is.na(matexSa))==0)) cat(sum(rowSums(is.na(matexSa))==0)," samples found everywhere\n")
  ##
  meta=.joinDF2(lapply(re,function(x) x$Meta),c("Sid","sType"),matexSa)
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
  
  ##############
  lso=order(meta$sType,meta$InjOrder)
  meta=meta[lso,]
  lusamp=lusamp[lso]
  fileinfos=fileinfos[lso,]
  matexSa=sapply(re,function(x) match(lusamp,x$Sid))
  
  ##############
  cat("Check analytes:\n")
  lumet=unlist(lapply(re,function(x) x$Analyte))
  lann=lapply(re,function(x) x$Annot)#[,grep("^[A-Za-z]+$",names(x$Annot))])
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
  cat("Merge data: ")
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
  cat("\n")
  
  allmat=lapply(allmat,function(x){dimnames(x)=list(meta$Sid,annot$Analyte);x})
  
  allmat=list(Method=names(re),Sid=meta$Sid,Analyte=annot$Analyte,Annot=annot,Meta=meta,File=fileinfos,Data=allmat)
  class(allmat)=append(class(allmat),"metaboSet")
  if(any(sapply(re,function(x) "fluxoSet"%in%class(x)))) class(allmat)=append(class(allmat),"fluxoSet")
  
  invisible(allmat)
}
