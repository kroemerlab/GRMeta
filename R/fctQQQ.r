
#  dictSyno=data.frame(Before=c("Dimethyl oxoglutaric acid","1-5.anhydroglucitol","Glyceraldehyde-3-phosphate" ,"Glyceraldehyde 3.phosphate", "Maltose Trehalose",
#                               "Dimethyl glycine","Di methyl glycine","Alpha tocopherol","Tocopherol","5-aminovaleric","Aminovaleric", "Phosphoenol pyruvate",
#                               "cytidine","adenosine","xanthosine","AMPc","AMP.10C13.5N15","COA", "ACETYLCOA POS", "ACETYLCOENZYME A.2C13", "malonyl coA"  ,  "SuccCoA"),
#                      After=c("Dimethyl-oxoglutaric acid","1,5-anhydroglucitol","Glyceraldehyde 3-phosphate", "Glyceraldehyde 3-phosphate" ,"Maltose;Trehalose" ,        
#                              "Dimethylglycine","Dimethylglycine","Alpha-tocopherol","Alpha-tocopherol","5-aminovaleric acid" , "5-aminovaleric acid","Phosphoenolpyruvic acid",
#                              "Cytidine","Adenosine","Xanthosine","cAMP","AMP_ISTD","CoA", "AcetylCoA", "AcetylCoA_ISTD", "MalonylCoA"  ,  "SuccinylCoA"  ),
#                      stringsAsFactors = FALSE)
#  save(file="dictSyno.rda",dictSyno)
#  dictRegExpr=data.frame(Before=c("__","_","^Beta\\.","^O\\.","^N\\.","^X([0-9]+)\\.","([A-Za-z])\\.([0-9])\\.([A-Za-z])","([A-Za-z])\\.([A-Za-z])","\\.$"," C13"),
#                         After=c(";"," ","Beta-","O-","N-","\\1-","\\1-\\2-\\3","\\1 \\2","","_ISTD"),stringsAsFactors = FALSE)
#  save(file="dictRegExpr.rda",dictRegExpr)
 
 

cleanMetaboNames<-function(metnam,RegExpr=dictRegExpr,Syno=dictSyno){

   
  vnam=metnam
  if(is.data.frame(RegExpr) & !is.null(RegExpr))  for(i in 1:nrow(RegExpr)) vnam=gsub(RegExpr$Before[i],RegExpr$After[i],vnam)
  if(is.data.frame(Syno) & !is.null(Syno)) for(i in 1:nrow(Syno)) vnam[which(vnam==Syno$Before[i])]=Syno$After[i]
  
  return(list(newnam=vnam,oldnams=metnam))
}


#################################################################################################

addSBRMatrix<-function(obj,thresh=1,SidBlank=NULL,what="Height",newname="SBR",fct=median){
  
  lbl=which(obj$Sid%in%SidBlank)
  if(!what%in%names(obj$Data)){cat(what," doesn't exist\n");return(obj)}
  mat=obj$Data[[what]]
  if(length(lbl)==0 | is.null(SidBlank)) norm=rep(thresh,length(obj$Analyte))
  if(length(lbl)>0) norm=exp(apply(log(mat[lbl,,drop=F]),2,fct,na.rm=T))
  norm[is.na(norm) | norm<thresh]=thresh
  nmat=sweep(mat,2,norm,"/")
  nmat[which(nmat<1)]=1
  obj$Data[[newname]]=nmat
  obj$Annot$SBR=norm
  return(obj)
}





