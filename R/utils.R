### check if x falls in the range or adjust to fit 
.GRchkrange<-function(x,rangeSc=c(1,Inf)){
  if(is.vector(x)){
    x=range(x,na.rm=T)
    return(c(max(floor(x[1]),rangeSc[1]),min(ceiling(x[2]),rangeSc[2])))
  }
  x=t(apply(x,1,range))
  x[,1]= floor(x[,1]) %>% {.[.<rangeSc[1]] = rangeSc[1]; .} 
  x[,2]= ceiling(x[,2]) %>% {.[.>rangeSc[2]] = rangeSc[2]; .} 
  x
}

# .GRchkrange2<-function(x,rangeSc=c(1,1,Inf,Inf)){
#   if(is.vector(x)){
#     x=c(floor(min(x,na.rm=T)),ceiling(max(x,na.rm=T)))
#     if(!(x[2]>rangeSc[2] | x[1]<rangeSc[3])) return(c(NA,NA))
#     ifelse(x[1])
#        
#     return(c(max(max(x[1],rangeSc[1]),rangeSc[2]),min(x[2]),rangeSc[2])))
#   }
#   x=t(apply(x,1,range))
#   x[,1]= floor(x[,1]) %>% {.[.<rangeSc[1]] = rangeSc[1]; .} 
#   x[,2]= ceiling(x[,2]) %>% {.[.>rangeSc[2]] = rangeSc[2]; .} 
#   x
# }

##### 
printElapse<-function(strtime,bef="",add="\n"){
  if(typeof(strtime)=="integer"){
    top=as.integer(as.POSIXct( Sys.time() ))-strtime
    units="secs"
  } else{
    top=Sys.time()-strtime
    units=attr(top,"units")
  }
  sprintf("%s%.2f %s%s",bef,top,units,add)
}

##################################################################################################
.GRgetAgilentXMLInfos<-function(cefi){
  
  ### Slow bit using library(XML)
  system.time(doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))[[2]])
  inst=xmlChildren(xmlChildren(doc)$instrument)
  dproc=xmlChildren(xmlChildren(doc)$dataProcessing)
  top=do.call("cbind",lapply(c(xmlToList(inst$detector),xmlToList(inst$source),xmlToList(inst$additional),
                               xmlToList(dproc$processingMethod)),function(x) x[c("name","value")]))
  df=top["value",]
  names(df)=top["name",]
  df=c(unlist(xmlToList(xmlChildren(xmlChildren(doc)$admin)$sourceFile)),df)
  df=data.frame(t(as.matrix(df,nrow=1)),stringsAsFactors = FALSE)
  dts=xmlAttrs(dproc$software)["completionTime"]
  dts = strsplit(gsub("\"", "", dts), "T")[[1]]
  df$completionTime=chron(dts[1], dts[2], format = c(dates = "Y-M-D", times = "h:m:s"))
  
  ### XCMS for polarity etc...
  xraw=xcms:::xcmsRaw(cefi)
  df$polarity=as.character(xraw@polarity[1])
  df$rtmin=round(min(xraw@scantime)/60,3)
  df$rtmax=round(max(xraw@scantime)/60,3)
  df$nscan=length(xraw@scantime)
  df$mzmin=round(min(xraw@mzrange),6)
  df$mzmax=round(max(xraw@mzrange),6)
  
  rownames(df)=NULL
  
  invisible(df)
}

##################################################################################################
.GRgetReAdWXMLInfos<-function(cefi){
  
  #########################
  ## get completion date
  cefiraw=c(gsub("\\.mzXML",".raw",cefi),gsub("\\.mzXML",".Raw",cefi),gsub("\\.mzXML",".RAW",cefi))
  cefiraw=cefiraw[file.exists(cefiraw)]
  dts=NA
  if(length(cefiraw)>0) dts=as.chron(file.info(cefiraw)[,"mtime"], format = c(dates = "Y-M-D", times = "h:m:s"))
  
  ########################
  ## get run infos
  tt=mzR:::openMSfile(cefi)
  infos=unlist(mzR:::runInfo(tt))
  df=data.frame(fileName=fileName(tt),completionTime=dts,polarity=NA,
                rtmin=round(infos["dStartTime"]/60,3),rtmax=round(infos["dEndTime"]/60,3),
                nscan=infos['scanCount'],
                mzmin=round(infos["lowMz"],6),mzmax=round(infos["highMz"],6))
  infos=unlist(mzR:::instrumentInfo(tt))
  for(i in names(infos)) df[,i]=infos[i]
  mzR:::close(tt)
  
  ########################
  ## polarity
  xraw = xcms:::xcmsRaw(cefi)
  df$polarity = as.character(xraw@polarity[1])
  rm(list="xraw")
  
  rownames(df)=NULL
  invisible(df)
}
