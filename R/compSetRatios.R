
# 
# load("/home/david/Metabo/Packaging/GRMeta/data/RatiosGR.rda")
# RatiosGR$`AcSpmd/Spmd`=list(c("N8-acetylspermidine","N1-acetylspermidine","N-acetylspermidine"),c("Spermidine"),8)
# RatiosGR$`AcSpr/Spr`=list(c("N-acetylspermine","N12-acetylspermine","N1-acetylspermine"),c('Spermine'),8)
# RatiosGR$`DiAcSpmd/AcSpmd`=list('N1,N8-diacetylspermidine',c("N8-acetylspermidine","N1-acetylspermidine","N-acetylspermidine"),8)
# RatiosGR$`DiAcSpr/AcSpr`=list('N1,N12-diacetylspermine',c("N-acetylspermine","N8-acetylspermine","N12-acetylspermine"),8)
# save(file="/home/david/Metabo/Packaging/GRMeta/data/RatiosGR.rda",RatiosGR)



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
    #if(all(unlist(lpairs[[ix]])%in%cmet$PutNam)){
    
    lnums=which(cmet$PutNam%in%lpairs[[ix]][[1]])
    n1=length(lnums)
    if(length(lnums)==0) next
    ldens=which(cmet$PutNam%in%lpairs[[ix]][[2]])
    if(length(ldens)==0) next
    lnums=rep(lnums,length(ldens))
    ldens=rep(ldens,each=n1)
    
    if(length(ldens)>1 & any(cmet$Ds[ldens]==cmet$Ds[lnums])){
      l2k=which(cmet$Ds[ldens]==cmet$Ds[lnums])
      ldens=ldens[l2k];lnums=lnums[l2k]
    }
    
    #   addm=cmat[,lnums,drop=F]-cmat[,ldens,drop=F]
    NNam=gsub("_0","",paste(ix,(1:length(lnums)-1),sep="_"))
    meth=cbind(cmet$Ds[lnums],cmet$Ds[ldens])
    meth=apply(meth,1,function(x) paste(unique(x),collapse="/"))
    levann=apply(cbind(cmet$lvA[lnums],cmet$lvA[ldens]),1,max)
    icmb=data.frame(Analyte=NNam,MetName=ix,AnalyteNum=cmet$Id[lnums],AnalyteDen=cmet$Id[ldens],LevelAnnot=levann,Method=meth,stringsAsFactors=F)
    icmb$Fam=lpairs[[ix]][[3]]
    cmb[[ix]]=icmb
    #}
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
  
  if(!conv2set)   return(invisible(list(Mat=Data,Stats=ratstats)))
  
  allmat=list(Method=obj$Method,Sid=obj$Sid,Analyte=ratstats$Analyte,Annot=ratstats,Meta=obj$Meta,File=obj$File,Data=Data)
  class(allmat)=append(class(allmat),"metaboSet")
  
  invisible(allmat)
  
  }
