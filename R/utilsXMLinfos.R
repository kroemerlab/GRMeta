##################################################################################################
.GRgetDirXMLInfos<-function(dir,type="Thermo",params = list(),short=T){
  
  pattern = ifelse(type=="Thermo",'\\.mzXML$','\\.mzdata.xml$')
  paramsvals <- GRMeta::paramsParsing()
  if (!missing(params)) 
    paramsvals[names(params)] <- params
  params = paramsvals
  
  dir=dir[file.exists(dir)]
  if(length(dir)==0) stop('Dir/file do not exist!')
  lfi=unique(c(dir[grep(pattern,dir)],list.files(dir,pattern = pattern,full.names = T)))
  cat("Found ",length(lfi),ifelse(type=="Thermo"," Thermo/ReAdW"," Agilent/MH")," xml files in ",dir,"\n",sep="")
  if(length(lfi)==0) return(NULL)
  adf=list()
  for(ifi in lfi){
    cat(".")
    re<-try(.GRgetXMLFileInfos(ifi,type),TRUE)
    if( "try-error"%in%class(re)){
      cat("Error for",ifi,"\n")
      next
    }
    adf[[ifi]]=re
  }
  lvar=unique(unlist(sapply(adf,names)))
  for(ifi in names(adf)){
    if(all(lvar%in%names(adf[[ifi]]))) next
    itoadd=lvar[!lvar%in%names(adf[[ifi]])]
    for(i in itoadd) adf[[ifi]][,i]=NA
  }
  adf<-do.call("rbind",adf)
  ##### check if date is.na due to NEG/POS split
  l=which(is.na(adf$completionTime) & grepl(paste0("_[NEGPOS]{3}",pattern),adf$fileName))
  if(length(l)) for(i in l){
  l2=which(adf$fileName==gsub("_[POSNEG]{3}\\.",".",adf$fileName[i]) & adf$dirName==adf$dirName[i])
  if(length(l2)==1) adf$completionTime[i]=adf$completionTime[l2] 
  }
  
  #####
  filenam=adf$fileName
  nams = gsub("([^.]+)\\.[[:alnum:]\\.]+$", "\\1", filenam)
  if (!is.null(params$NameClean)) for (i in params$NameClean) nams = gsub(i, "", nams)
  styp = rep(NA, length(nams))
  if (length(grep(params$regTypes, nams)) > 0) 
    styp[grep(params$regTypes, nams)] = gsub(params$regTypes, "\\1", nams)[grep(params$regTypes, nams)]
  if (any(is.na(styp))) {
    cat(sum(is.na(styp)), "unknown sample type:\n", nams[which(is.na(styp))], 
        "\n", sep = " ")
    styp[is.na(styp)] = "Unknown"
  }
  sid = nams
  adf = cbind(data.frame(Sid = sid, sType = styp, InjOrder = order(order(adf$completionTime)), stringsAsFactors = FALSE),adf)
  l2exp=c("Sid","sType","InjOrder","dirName","fileName","completionTime",'acquisitionMethod',
       'msLevels', "rtmin","rtmax","polarity","nscan","mzmin","mzmax","polarity2",'nscan2',"mzmin2","mzmax2",'collisionEnergy')
  if(short) adf=adf[,names(adf)%in%l2exp]
  if(max(table(adf$Sid)==1)) rownames(adf)=adf$Sid
  return(adf)
}


##################################################################################################
.GRgetXMLFileInfos<-function(cefi,Type=ifelse(grepl("mzdata",cefi),"Agilent","Thermo")){
  
  require(mzR)
  require(xcms)
  
  ########################
  ## get run infos
  hasMS2=FALSE
  tt=mzR:::openMSfile(cefi)
  infos=unlist(mzR:::runInfo(tt))
  lmslevels=infos[grep('msLevels',names(infos))]
  ttheader=header(tt)
  if(any(lmslevels==2)) hasMS2=T
  df=data.frame(dirName=normalizePath(dirname(cefi)),
                fileName=basename(cefi),completionTime=NA,polarity=NA,msLevels=paste(lmslevels,collapse="/"),
                rtmin=round(infos["dStartTime"]/60,4),rtmax=round(infos["dEndTime"]/60,4),
                nscan=sum(ttheader$msLevel==1),
                mzmin=round(min(ttheader$lowMZ[ttheader$msLevel==1]),6),
                mzmax=round(max(ttheader$highMZ[ttheader$msLevel==1],6)))
  if(any(lmslevels==2)){
    hasMS2=T
    df$nscan2=sum(ttheader$msLevel==2)
    df$mzmin2=round(min(ttheader$lowMZ[ttheader$msLevel==2]),6)
    df$mzmax2=round(max(ttheader$highMZ[ttheader$msLevel==2]),6)
  }
  infos=unlist(mzR:::instrumentInfo(tt))
  infos=infos[infos!=""]
  for(i in names(infos)) df[,i]=infos[i]
  mzR:::close(tt)
  
  ########################
  ## polarity
  suppressWarnings(xraw<-xcms:::xcmsRaw(cefi))
  df$polarity=ifelse(length(xraw@polarity)==0,"none",paste(unique(as.character(xraw@polarity)),collapse="/"))
  rm(list="xraw")
  if(hasMS2){
    suppressWarnings(xraw<-xcms:::xcmsRaw(cefi,includeMSn = T))
    df$polarity2=ifelse(length(xraw@polarity)==0,"none",paste(unique(as.character(xraw@polarity)),collapse="/"))
    if(Type=="Thermo"){
      doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))
      allatr=do.call("rbind",lapply(xmlChildren(doc[[1]]),function(x) do.call("rbind",lapply(xmlChildren(x),function(y) xmlAttrs(y)[c("msLevel","filterLine")]))))
      allatr=allatr[which(allatr[,1]=="2"),]
      lce=paste(unique(gsub(".*[0-9]+@(.*) .*","\\1",allatr[,2])),collapse="/")
      df$collisionEnergy=lce
      rm(list="doc")
      
    }
  }
  
  ###################################
  ## get completion date for Thermo
  if(Type=="Thermo"){
    cefiraw=paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", cefi),c(".raw",".Raw",".RAW"))
    cefiraw=cefiraw[file.exists(cefiraw)]
    if(length(cefiraw)>0) df$completionTime=as.chron(file.info(cefiraw)[,"mtime"], format = c(dates = "Y-M-D", times = "h:m:s"))
  }
  if(Type=="Agilent"){
    ### Slow bit using library(XML)
    require(XML)
    system.time(doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))[[2]])
    inst=xmlChildren(xmlChildren(doc)$instrument)
    dproc=xmlChildren(xmlChildren(doc)$dataProcessing)
    top=do.call("cbind",lapply(c(xmlToList(inst$detector),xmlToList(inst$source),xmlToList(inst$additional),
                                 xmlToList(dproc$processingMethod)),function(x) x[c("name","value")]))
    toadd=top["value",]
    names(toadd)=top["name",]
    toadd=c(unlist(xmlToList(xmlChildren(xmlChildren(doc)$admin)$sourceFile)),toadd)
    toadd=data.frame(t(as.matrix(toadd,nrow=1)),stringsAsFactors = FALSE)
    dts=xmlAttrs(dproc$software)["completionTime"]
    dts = strsplit(gsub("\"", "", dts), "T")[[1]]
    df$completionTime=chron(dts[1], dts[2], format = c(dates = "Y-M-D", times = "h:m:s"))
    df=cbind(df,toadd)
  }
  
  rownames(df)=NULL
  return(df)
}


# ##################################################################################################
# .GRgetAgilentXMLInfos<-function(cefi){
#   
#   ### Slow bit using library(XML)
#   system.time(doc <- xmlRoot(xmlParse(cefi, useInternal = TRUE))[[2]])
#   inst=xmlChildren(xmlChildren(doc)$instrument)
#   dproc=xmlChildren(xmlChildren(doc)$dataProcessing)
#   top=do.call("cbind",lapply(c(xmlToList(inst$detector),xmlToList(inst$source),xmlToList(inst$additional),
#                                xmlToList(dproc$processingMethod)),function(x) x[c("name","value")]))
#   df=top["value",]
#   names(df)=top["name",]
#   df=c(unlist(xmlToList(xmlChildren(xmlChildren(doc)$admin)$sourceFile)),df)
#   df=data.frame(t(as.matrix(df,nrow=1)),stringsAsFactors = FALSE)
#   dts=xmlAttrs(dproc$software)["completionTime"]
#   dts = strsplit(gsub("\"", "", dts), "T")[[1]]
#   df$completionTime=chron(dts[1], dts[2], format = c(dates = "Y-M-D", times = "h:m:s"))
#   
#   ### XCMS for polarity etc...
#   xraw=xcms:::xcmsRaw(cefi)
#   df$polarity=as.character(xraw@polarity[1])
#   df$rtmin=round(min(xraw@scantime)/60,3)
#   df$rtmax=round(max(xraw@scantime)/60,3)
#   df$nscan=length(xraw@scantime)
#   df$mzmin=round(min(xraw@mzrange),6)
#   df$mzmax=round(max(xraw@mzrange),6)
#   
#   df=cbind(data.frame(dirName=normalizePath(dirname(cefi)),
#                       fileName=basename(cefi),stringsAsFactors = F),df)
#   rownames(df)=NULL
#   
#   invisible(df)
# }
