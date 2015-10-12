exportDataXLSX<-function(obj,outfile,ldata=names(obj$Data),sortSample=c("sType","InjOrder"),sortAnalyte=c("Method","RT"),transpose=FALSE,nround=3,characterNA=c("-","NA")){
  
  if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
  
  ldata=ldata[which(ldata%in%names(obj$Data))]
  sortAnalyte=sortAnalyte[which(sortAnalyte%in%names(obj$Annot))]
  sortSample=sortSample[which(sortSample%in%names(obj$Meta))]
  
  lsoAna=1:length(obj$Analyte)
  if(!is.null(sortAnalyte)){
    lsoAna=paste("order(",paste(paste("obj$Annot$",sortAnalyte,sep=""),collapse=","),")",sep="")
    lsoAna=eval(parse(text=lsoAna))
  }
  
  lsoSa=1:length(obj$Sid)
  if(!is.null(sortSample)){
    lsoSa=paste("order(",paste(paste("obj$Meta$",sortSample,sep=""),collapse=","),")",sep="")
    lsoSa=eval(parse(text=lsoSa))
  }
  
  indf<-function(df,nround,characterNA,what="Var"){
    for(j in names(df)){
      if(is.numeric(df[,j])) df[,j]=round(df[,j],nround)
      df[,j]=as.character(df[,j])
      if(any(is.na(df[,j]))) df[is.na(df[,j]),j]=characterNA 
    }
   df=cbind("Var0000000"=as.character(rownames(df)),df,stringsAsFactors=FALSE)
   names(df)[1]=what
#   print(names(df))
    df=rbind(names(df),df)
    df
  }
  
  cat("Exporting to ",outfile," :\n",sep="")
  write.xlsx2(indf(obj$Meta[lsoSa,which(names(obj$Meta)!="Sid")],nround,characterNA[1],"Sid"),
                file = outfile,col.names=FALSE, row.names=FALSE,append=FALSE,sheetName ="SampleInfos")  
  write.xlsx2(indf(obj$File[lsoSa,which(names(obj$File)!="Sid")],nround,characterNA[1],"Sid"),
                file = outfile,col.names=FALSE, row.names=FALSE,append=TRUE,sheetName ="FileInfos")  
  write.xlsx2(indf(obj$Annot[lsoAna,which(names(obj$Annot)!="Analyte")],nround,characterNA[1],"Analyte"),
                file = outfile,col.names=FALSE, row.names=FALSE,append=TRUE,sheetName ="VarInfos")
  
  for(i in ldata){
    if(transpose) m=round(t(obj$Data[[i]][lsoSa,lsoAna,drop=FALSE]),nround)
    if(!transpose) m=round(obj$Data[[i]][lsoSa,lsoAna,drop=FALSE],nround)
    m=apply(m,2,as.character)
#    print(str(m))
    if(transpose) df=cbind("Analyte"=obj$Analyte[lsoAna],m)
    if(!transpose) df=cbind("Sid"=obj$Sid[lsoSa],m)
    #   print(names(df))
    df=rbind(colnames(df),df)
    df[is.na(df)]=characterNA[2]
    
    write.xlsx2(data.frame(df,stringsAsFactors = FALSE),
                    file = outfile,col.names=FALSE, row.names=FALSE,append=TRUE,sheetName =i)
    cat(i," ",sep="")
  }
  
  invisible(list(lsoSa,lsoAna))
}
