
compSetRatios<-function(obj,l2add=names(obj$Data),lpairs=RatiosGR,conv2set=FALSE){
  
#   if(is.null(lpairs)){
#     RatiosGR<-NULL
#     data("RatiosGR")
#     lpairs=RatiosGR
#     rm("RatiosGR")
#   }
  
  lnna=which(!is.na(obj$Annot$MetName))
  cmet=do.call("rbind",lapply(lnna,function(ix) data.frame(Id=obj$Annot$Analyte[ix],Typ=1,Ds=obj$Annot$Method[ix],lvA=obj$Annot$LevelAnnot[ix],
                                                           PutNam=strsplit(obj$Annot$MetName[ix],"[;/]")[[1]],stringsAsFactors=FALSE)))
  cmet$PutNam=gsub("\\?","",cmet$PutNam)
  
  
  cmb=list()
  for(ix in names(lpairs)){
    if(all(unlist(lpairs[[ix]])%in%cmet$PutNam)){
      
      lnums=which(cmet$PutNam==lpairs[[ix]][[1]]);n1=length(lnums)
      ldens=which(cmet$PutNam==lpairs[[ix]][[2]])
      lnums=rep(lnums,length(ldens))
      ldens=rep(ldens,each=n1)
      
      if(length(ldens)>1 & any(cmet$Ds[ldens]==cmet$Ds[lnums])){
        l2k=which(cmet$Ds[ldens]==cmet$Ds[lnums])
        ldens=ldens[l2k];lnums=lnums[l2k]
      }
      
      #   addm=cmat[,lnums,drop=F]-cmat[,ldens,drop=F]
      NNam=gsub("_0","",paste(ix,(1:length(lnums)-1),sep="_"))
      meth=cbind(cmet$Ds[lnums],cmet$Ds[ldens])
      meth=apply(meth,1,function(x) paste(unique(sort(x)),collapse="/"))
      levann=apply(cbind(cmet$lvA[lnums],cmet$lvA[ldens]),1,max)
      icmb=data.frame(Analyte=NNam,MetName=ix,AnalyteNum=cmet$Id[lnums],AnalyteDen=cmet$Id[ldens],LevelAnnot=levann,Method=meth,stringsAsFactors=F)
      
      cmb[[ix]]=icmb
    }
  }
  ratstats=do.call("rbind",cmb)
  rownames(ratstats)=ratstats$Analyte
  
  Data=list()
  for(i in l2add[l2add%in%names(obj$Data)]){
    mat=obj$Data[[i]]
    if(grepl("^Log",i))
      Data[[i]]=do.call("cbind",lapply(rownames(ratstats),function(x) mat[,ratstats[x,]$AnalyteNum]-mat[,ratstats[x,]$AnalyteDen]))
    if(!grepl("^Log",i))
      Data[[i]]=do.call("cbind",lapply(rownames(ratstats),function(x) mat[,ratstats[x,]$AnalyteNum]/mat[,ratstats[x,]$AnalyteDen]))
    colnames(Data[[i]])=rownames(ratstats)
  }
  
  if(conv2set){
    allmat=list(Method=obj$Method,Sid=obj$Sid,Analyte=ratstats$Analyte,Annot=ratstats,Meta=obj$Meta,File=obj$File,Data=Data)
    class(allmat)=append(class(allmat),"metaboSet")
  } else {allmat=list(Mat=Data,Stats=ratstats)}
  
  invisible(allmat)
  
}
