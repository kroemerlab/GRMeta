print.metaboSet<-function(x,length=12){
  
  cat("** ",length(x$Analyte)," analyte",ifelse(length(x$Analyte)==1,"","s"),": ",paste(names(table(x$Annot$Method)),"(",(table(x$Annot$Method)),") ",sep=""),"\n",sep="")
  cat("** Annotation level: ",paste(names(table(x$Annot$LevelAnnot)),"(",(table(x$Annot$LevelAnnot)),") ",sep=""),"\n",sep="")
  n=min(length(x$Analyte),length)
  top=x$Analyte[1:n]
  add=length(x$Analyte)-n;if(add>0) top=c(top,paste("+",add,"entries"))
#  print.default(top)
#cat("\n")
  cat("** ",length(x$Sid)," sample",ifelse(length(x$Sid)==1,"","s"),": ",paste(names(table(x$Meta$sType)),"(",(table(x$Meta$sType)),") ",sep=""),"\n",sep="")
  n=min(length(x$Sid),length)
  top=x$Sid[1:n]
  add=length(x$Sid)-n;if(add>0) top=c(top,paste("+",add,"entries"))
#  print.default(top)
#  cat("\n")
  cat("** ",length(x$Data)," data matri",ifelse(length(x$Data)==1,"x","ces"),": ",paste(names(x$Data),collapse=" "),"\n",sep="")
  
}

