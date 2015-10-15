compSetSBR<-function(obj,thresh=1,SidBlank=NULL,what="Height",newname=NULL,fct=median){
  
  lbl=which(obj$Sid%in%SidBlank)
  if(!what%in%names(obj$Data)){cat(what," doesn't exist\n");return(obj)}
  mat=obj$Data[[what]]
  if(length(lbl)==0 | is.null(SidBlank)) norm=rep(thresh,length(obj$Analyte))
  if(length(lbl)>0) norm=exp(apply(log(mat[lbl,,drop=F]),2,fct,na.rm=T))
  norm[is.na(norm) | norm<thresh]=thresh
  nmat=sweep(mat,2,norm,"/")
  nmat[which(nmat<1)]=1
  
  if(is.null(newname)) invisible(list(SBR=nmat,Norm=norm)) else{
    obj$Data[[newname]]=nmat
    obj$Annot$SBR=norm
    invisible(obj)
  }
  
}
