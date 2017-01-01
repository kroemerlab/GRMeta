
.GRparseMRMMethod<-function(methfile){
  dat <- xmlRoot(xmlTreeParse(methfile))
  tSegs=dat[['timeSegments']][[1]][['scanSegments']]
  SegRes=list()
  for(x in xmlChildren(tSegs)){
    idx=xmlValue(x[["index"]])
    SegRes[[idx]]=lapply(xmlChildren(x[["scanElements"]]),function(y)
      c(SegIdx=idx,SegIdx2=xmlValue(y[["index"]]),Q1=xmlValue(y[["ms1LowMz"]]),Q3=xmlValue(y[["ms2LowMz"]]),Cmpd=xmlValue(y[["compoundName"]]),
        Frag=xmlValue(y[["fragmentor"]]),Dwell=xmlValue(y[["dwell"]]),CE=xmlValue(y[["collisionEnergy"]]),EMV=xmlValue(y[["deltaEMV"]])))
  }
  suppressWarnings(SegDF<-data.frame(do.call("rbind",unlist(SegRes,recursive = F))))
  for(i in c("SegIdx","SegIdx2","Q1","Q3","Dwell","Frag","CE","EMV")) SegDF[,i]=as.numeric(SegDF[,i])
  SegDF$MRM=paste0(SegDF$Q1,">",SegDF$Q3)
  SegDF$MRMId=paste0(SegDF$MRM,"_",SegDF$Frag,"_",SegDF$CE)
  if(any(SegDF$SegIdx2>1)) SegDF$MRMId=paste0(SegDF$MRM,"_",SegDF$Frag,"_",SegDF$CE,"-",SegDF$SegIdx2)
  rownames(SegDF)=SegDF$MRMId
  invisible(SegDF)
}

##### xmlfile

parseAgilentMRMXML<-function(xmlfile,methinfos,doDate=TRUE,commonRT=NULL){
  
  #### 
  if(!is.data.frame(methinfos)){
    if(!file.exists(methinfos)) stop(paste0('Not a dataframe and the method file does not exists!!!'))
    cat('Parsing the method file ',methinfos,"\n",sep="")
    methinfos=.GRparseMRMMethod(methinfos)
  }
  #### 
  cat('Parsing the data file ',xmlfile,"\n",sep="")
  cdat=mzR:::openMSfile(xmlfile)
  hcdat=mzR:::header(cdat)
  hcdat=hcdat[,!apply(hcdat,2,function(x) all(x==0))]
  msdat=mzR::peaks(cdat)
  msdat=do.call("rbind",lapply(1:length(msdat),function(i) cbind(sc=i,msdat[[i]])))
  msdat=cbind(msdat,rt=hcdat[msdat[,1],"retentionTime"]/60)
  colnames(msdat)[2:3]=c("mz","y")
  
  n=floor(nrow(msdat)/nrow(methinfos))
  msdat=msdat[1:(n*nrow(methinfos)),]
  hcdat=hcdat[hcdat[,1]%in%msdat[,1],]
  hcdat=hcdat[match(msdat[,"sc"],hcdat[,1]),]
  lso=order(methinfos$SegIdx,methinfos$Q3)
  msdat=cbind(mrm=rep(lso,n),scid=rep(1:n,each=length(lso)),msdat)
  
  tt=table(hcdat[,1],hcdat[,'peaksCount'])
  chk1=all(sapply(colnames(tt),function(i) all(tt[,i]%in%c(0,as.numeric(i)))))
  chk2=all(tapply(msdat[,"mz"],msdat[,1],function(x) length(unique))==1)
  chk3=all(methinfos[msdat[,"mrm"],"Q3"]==msdat[,"mz"])
  
  if(!(chk1 & chk2 & chk3)){
    cat(" ********* WARNINGS *********\n")
    cat(" --- peaksCount and number in spectra match:",chk1,"\n")
    cat(" --- spectra mz and MRM id match:",chk2,"\n")
    cat(" --- spectra mz and MRM Q3 match:",chk3,"\n")
    return(NULL)
  }
  
  cat(" -->  it looks OK!\n")
  cat(" *** MRM: ",nrow(methinfos),"\n",sep="")
  cat("    - unique MRM=",length(unique(methinfos$MRM))," / cmpds=",length(unique(methinfos$Cmpd)),"\n",sep="")
  cat("    - unique Q1=",length(unique(methinfos$Q1))," / Q3=",length(unique(methinfos$Q3)),"\n",sep="")
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
  dimnames(aint)=dimnames(art)=list(luscid,rownames(methinfos)[as.numeric(lumrm)])
  
  cDate=ifelse(doDate,.GRgetcomptime(xmlfile),NA)
  
  if(!is.null(commonRT)) if(any(is.na(commonRT))) return(list(MSdat=msdat,RT=art,Int=aint,IntRT=NULL,SegDF=methinfos,Date=cDate))
  
  ## rough estimate of a common retention time
  if(is.null(commonRT)){
  rrt=c(max(apply(art,2,min)),min(apply(art,2,max)))
  ## should be quite similar
  table(colSums(art>rrt[1] & art<rrt[2]))
  ## sexier drt???
  drt=apply(art,2,function(x) median(diff(x)))
  drt=median(round(drt,ceiling(-log10(drt))+2))
  drt=round(drt,ceiling(-log10(drt))+2)
  ## set new rt scale
  newdrt=c(max(ceiling(rrt[1]/drt),1),floor(rrt[2]/drt)-1)*drt
  table(colSums(art>=newdrt[1] & art<=newdrt[2]))
  commonRT=c(newdrt,drt)
  cat(sprintf(" *** common estimated RT from %.4f to %.4f by %.4f min.\n",commonRT[1],commonRT[2],commonRT[3]))
  }
  if(length(commonRT)>3)  newrt=commonRT else newrt=seq(commonRT[1],commonRT[2],commonRT[3])
 
  ## linear approximation -> cubic spline instead
  aint2=lapply(1:ncol(aint),function(ix) approx(x=art[,ix],y=aint[,ix],xout = newrt)$y)
  aint2=do.call("cbind",aint2)
  dimnames(aint2)=list(newrt,colnames(aint))
 
  invisible(list(MSdat=msdat,RT=art,Int=aint,IntRT=aint2,SegDF=methinfos,Date=cDate))
  
  
}

