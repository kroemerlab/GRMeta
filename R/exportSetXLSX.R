
exportSetXLSX<-function(obj,outfile,ldata=names(obj$Data),sortSample=c("sType","InjOrder"),sortAnalyte=c("Method","RT"),transpose=FALSE,nround=3,characterNA=c("-","NA")){
  
  jgc <- function(){
    .jcall("java/lang/System", method = "gc")
  }   
  
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
  
  wb <- createWorkbook()  
  cat("Exporting to ",outfile," :\n",sep="")
  toadd=indf(obj$Meta[lsoSa,which(names(obj$Meta)!="Sid")],nround,characterNA[1],"Sid")
  jgc()
  cat("SampleInfos ")
  sheet <- createSheet(wb, sheetName = "SampleInfos")
  addDataFrame(toadd,sheet, col.names=FALSE, row.names=FALSE)

  toadd=indf(obj$File[lsoSa,which(names(obj$File)!="Sid")],nround,characterNA[1],"Sid")
  jgc()
  cat("FileInfos ")
  sheet <- createSheet(wb, sheetName = "FileInfos")
  addDataFrame(toadd,sheet, col.names=FALSE, row.names=FALSE)

  toadd=indf(obj$Annot[lsoAna,which(names(obj$Annot)!="Analyte")],nround,characterNA[1],"Analyte")
  jgc()
  cat("VarInfos ")
  sheet <- createSheet(wb, sheetName = "VarInfos")
  addDataFrame(toadd, sheet, col.names=FALSE, row.names=FALSE)
  
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
    
    cat(i," ",sep="")
    jgc()
    sheet <- createSheet(wb, sheetName = i)
    addDataFrame(df,sheet, col.names=FALSE, row.names=FALSE)
  
   }
  saveWorkbook(wb, outfile) 
  invisible(list(lsoSa,lsoAna))
}
# 
# exportSetXLSX2<-function(obj,outfile,ldata=names(obj$Data),sortSample=c("sType","InjOrder"),sortAnalyte=c("Method","RT"),transpose=FALSE,nround=3,characterNA=c("-","NA")){
#   
# # ldata=names(obj$Data);sortSample=c("sType","InjOrder");sortAnalyte=c("Method","RT");transpose=FALSE;nround=3;characterNA=c("-","NA")
#   
#   if(!inherits(obj, "metaboSet")) stop("This is not a metaboSet object")
#   
#   ldata=ldata[which(ldata%in%names(obj$Data))]
#   sortAnalyte=sortAnalyte[which(sortAnalyte%in%names(obj$Annot))]
#   sortSample=sortSample[which(sortSample%in%names(obj$Meta))]
#   
#   lsoAna=1:length(obj$Analyte)
#   if(!is.null(sortAnalyte)){
#     lsoAna=paste("order(",paste(paste("obj$Annot$",sortAnalyte,sep=""),collapse=","),")",sep="")
#     lsoAna=eval(parse(text=lsoAna))
#   }
#   
#   lsoSa=1:length(obj$Sid)
#   if(!is.null(sortSample)){
#     lsoSa=paste("order(",paste(paste("obj$Meta$",sortSample,sep=""),collapse=","),")",sep="")
#     lsoSa=eval(parse(text=lsoSa))
#   }
#   
#   indf<-function(df,nround,characterNA,what="Var"){
#     for(j in names(df)){
#       if(is.numeric(df[,j])) df[,j]=round(df[,j],nround)
#       df[,j]=as.character(df[,j])
#       if(any(is.na(df[,j]))) df[is.na(df[,j]),j]=characterNA 
#     }
#     df=cbind("Var0000000"=as.character(rownames(df)),df,stringsAsFactors=FALSE)
#     names(df)[1]=what
#     #   print(names(df))
#     df=rbind(names(df),df)
#     df
#   }
#   
#   cat("Exporting to ",outfile," :\n",sep="")
#   df1=indf(obj$Meta[lsoSa,which(names(obj$Meta)!="Sid")],nround,characterNA[1],"Sid");shn ="SampleInfos"
#   df2=indf(obj$File[lsoSa,which(names(obj$File)!="Sid")],nround,characterNA[1],"Sid");shn=append(shn,"FileInfos")
#   df3=indf(obj$Annot[lsoAna,which(names(obj$Annot)!="Analyte")],nround,characterNA[1],"Analyte");shn=append(shn,"VarInfos")
#   
#   j=3
#   for(i in ldata){
#     j=j+1
#     if(transpose) m=round(t(obj$Data[[i]][lsoSa,lsoAna,drop=FALSE]),nround)
#     if(!transpose) m=round(obj$Data[[i]][lsoSa,lsoAna,drop=FALSE],nround)
#     m=apply(m,2,as.character)
#     #    print(str(m))
#     if(transpose) df=cbind("Analyte"=obj$Analyte[lsoAna],m)
#     if(!transpose) df=cbind("Sid"=obj$Sid[lsoSa],m)
#     #   print(names(df))
#     df=rbind(colnames(df),df)
#     df[is.na(df)]=characterNA[2]
#     eval(expr = parse(text=paste("df",j,"=df",sep="")))
#     shn=append(shn,i)
# #     write.xlsx2(data.frame(df,stringsAsFactors = FALSE),
# #                 file = outfile,col.names=FALSE, row.names=FALSE,append=TRUE,sheetName =i)
#     cat(i," ",sep="")
#   }
#   shn=paste(shn,collapse = ',')
#   ldf=paste(paste("df",1:j,sep=""),collapse=",")
#   eval(expr = parse(text=paste("dataframes2xls:::write.xls(c(",ldf,"), outfile, sh.names =shn,col.names=FALSE,row.names=FALSE)",sep="")))
#   
#   invisible(list(lsoSa,lsoAna))
# }
