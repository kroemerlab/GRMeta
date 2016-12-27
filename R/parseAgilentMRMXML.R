
.GRparseMRMMethod<-function(methfile){
  dat <- xmlRoot(xmlTreeParse(methfile))
  tSegs=dat[['timeSegments']][[1]][['scanSegments']]
  SegRes=list()
  for(x in xmlChildren(tSegs)){
    idx=xmlValue(x[["index"]])
    SegRes[[idx]]=lapply(xmlChildren(x[["scanElements"]]),function(y)
      c(SegIdx=idx,SegIdx2=xmlValue(y[["index"]]),Q1=xmlValue(y[["ms1LowMz"]]),Q3=xmlValue(y[["ms2LowMz"]]),Cmpd=xmlValue(y[["compoundName"]]),
        Frag=xmlValue(y[["fragmentor"]]),Dwell=xmlValue(y[["dwell"]]),CE=xmlValue(y[["collisionEnergy"]])))
  }
  suppressWarnings(SegDF<-data.frame(do.call("rbind",unlist(SegRes,recursive = F))))
  for(i in c("SegIdx","SegIdx2","Q1","Q3","Dwell","Frag")) SegDF[,i]=as.numeric(SegDF[,i])
  SegDF$MRM=paste0(SegDF$Q1,">",SegDF$Q3)
  rownames(SegDF)=SegDF$MRMId=paste0(SegDF$MRM,"_",SegDF$Frag)
  invisible(SegDF)
}

parseAgilentMRMXML<-function(xmlfile,segdf,doDate=TRUE){
  
  #### 
  if(!is.data.frame(segdf)){
    if(!file.exists(segdf)) stop(paste0('Not a dataframe and the method file does not exists!!!'))
    cat('Parsing the method file ',segdf,"\n",sep="")
    segdf=.GRparseMRMMethod(segdf)
  }
  #### 
  cdat=mzR:::openMSfile(xmlfile)
  hcdat=mzR:::header(cdat)
  hcdat=hcdat[,!apply(hcdat,2,function(x) all(x==0))]
  msdat=mzR::peaks(cdat)
  msdat=do.call("rbind",lapply(1:length(msdat),function(i) cbind(sc=i,msdat[[i]])))
  msdat=cbind(msdat,rt=hcdat[msdat[,1],"retentionTime"]/60)
  colnames(msdat)[2:3]=c("mz","y")
  
  n=floor(nrow(msdat)/nrow(segdf))
  msdat=msdat[1:(n*nrow(segdf)),]
  hcdat=hcdat[hcdat[,1]%in%msdat[,1],]
  hcdat=hcdat[match(msdat[,"sc"],hcdat[,1]),]
  lso=order(segdf$SegIdx,segdf$Q3)
  msdat=cbind(mrm=rep(lso,n),scid=rep(1:n,each=length(lso)),msdat)
  
  tt=table(hcdat[,1],hcdat[,'peaksCount'])
  chk1=all(sapply(colnames(tt),function(i) all(tt[,i]%in%c(0,as.numeric(i)))))
  chk2=all(tapply(msdat[,"mz"],msdat[,1],function(x) length(unique))==1)
  chk3=all(segdf[msdat[,"mrm"],"Q3"]==msdat[,"mz"])
  
  if(!(chk1 & chk2 & chk3)){
    cat(" ********* WARNINGS *********\n")
    cat(" --- peaksCount and number in spectra match:",chk1,"\n")
    cat(" --- spectra mz and MRM id match:",chk2,"\n")
    cat(" --- spectra mz and MRM Q3 match:",chk3,"\n")
    return(NULL)
  }
  
  cat(xmlfile," looks OK:\n")
  cat(" *** MRM: ",nrow(segdf),"\n",sep="")
  cat("    - unique MRM=",length(unique(segdf$MRM))," / cmpds=",length(unique(segdf$Cmpd)),"\n",sep="")
  cat("    - unique Q1=",length(unique(segdf$Q1))," / Q3=",length(unique(segdf$Q3)),"\n",sep="")
  cat(" *** RT: ",sprintf("%.4f - %.4f",min(msdat[,"rt"]),max(msdat[,"rt"])))
  cat(" / num scan per MRM=",max(msdat[,"scid"]),"\n",sep="")
  
  art<-aint<-list()
  luscid=unique(msdat[,"scid"])
  lumrm=sort(unique(msdat[,"mrm"]))
  for(j in lumrm){
    l=which(msdat[,"mrm"]==j)
    l=l[match(luscid,msdat[l,"scid"])]
    art[[j]]=msdat[l,"rt"]
    aint[[j]]=msdat[l,"y"]
  }
  aint=do.call("cbind",aint)
  art=do.call("cbind",art)
  dimnames(aint)=dimnames(art)=list(luscid,rownames(segdf)[as.numeric(lumrm)])
  
  Date=ifelse(doDate,.GRgetcomptime(xmlfile),NULL)
  
  
  invisible(list(MSdat=msdat,RT=art,Int=aint,SegDF=segdf,Date=Date))
  
  
}

