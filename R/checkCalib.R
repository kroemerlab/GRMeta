.removeDups<-function(mdiffcal){
  l=which(mdiffcal[,"dups"]==1)
  l2chk=sort(unique(mdiffcal[l,"cal"]))
  for(i2chk in l2chk){
    lndups=which(mdiffcal[,"cal"]==i2chk & mdiffcal[,"dups"]==0)
    ldups=which(mdiffcal[,"cal"]==i2chk & mdiffcal[,"dups"]==1)
    cldups=.GRsplist(mdiffcal[ldups,1],ldups,d=2.1)
    lbest=unlist(lapply(cldups,function(y){
      ir=range(mdiffcal[y,1])
      lbef=lndups[which((ir[1]-mdiffcal[lndups,1])>0)];if(length(lbef)) lbef=lbef[1:min(5,length(lbef))]
      laft=lndups[which((mdiffcal[lndups,1]-ir[2])>0)];if(length(laft)) laft=laft[1:min(5,length(laft))]
      medaround=median(c(mdiffcal[lbef,"dppm"],mdiffcal[laft,"dppm"]))
      tapply(y,mdiffcal[y,"scan"],function(x) x[which.min(abs(mdiffcal[x,"dppm"]-medaround))])
    }))
    mdiffcal[lbest,"dups"]=0
  }
  mdiffcal}


chkCalib<-function(ifile,lcalib,maxdppm=121,maxdmz=0.001,minpts=5,rtlim=c(NA,NA)){
  
  xRaw <- xcmsRaw(ifile)
  
  ##### get all ions within maxdppm 
  mdiffcal=do.call("rbind",lapply(1:length(lcalib),function(ical){
    cat(lcalib[ical],"/",sep="")
    dmz=max(c(lcalib[ical]*maxdppm*10^-6,maxdmz),na.rm=T)
    x=.GRrawMat(xRaw,mzrange=lcalib[ical]*(1+1.001*c(-1,1)*maxdppm*10^-6))
    if(nrow(x)<minpts) return(NULL)
    isdups=x[,1]%in%names(which(table(x[,1])>1))*1
    cbind(x,cal=ical,calmz=lcalib[ical],dups=isdups)}))
  ## computing stats
  mdiffcal=cbind(mdiffcal,"dmz"=mdiffcal[,"mz"]-mdiffcal[,"calmz"])
  mdiffcal=cbind(mdiffcal,"dppm"=10^6*mdiffcal[,"dmz"]/mdiffcal[,"calmz"])
  mdiffcal0=mdiffcal=mdiffcal[order(mdiffcal[,1],mdiffcal[,"calmz"],abs(mdiffcal[,"dmz"])),]
  
  #######
  if(!is.na(rtlim[1])) mdiffcal=mdiffcal[mdiffcal[,"rt"]>=rtlim[1],,drop=F]
  if(!is.na(rtlim[2])) mdiffcal=mdiffcal[mdiffcal[,"rt"]<=rtlim[2],,drop=F]
  
  cat("\nFound ",length(unique(mdiffcal[,"cal"])),"/",length(lcalib)," calib m/z in ",length(unique(mdiffcal[,1])),"/",length(xRaw@scantime)," scans:\n",sep="")
  l=which(mdiffcal[,"dups"]==1)
  if(length(l)>0){
    x=rev(sort(rowSums(table(mdiffcal[l,"calmz"],mdiffcal[l,"scan"])>0)))
    cat("   ---> scans with duplicated calib: ",sum(x)," = ",paste0(round(as.numeric(names(x)),4),"(",x,")"),"\n")
    mdiffcal=.removeDups(mdiffcal0)
    mdiffcal=mdiffcal[which(mdiffcal[,"dups"]!=1),]
  }
  x=as.vector(table(mdiffcal[,1]))
  cat("   ---> num. calib per scan: ",sprintf("mean=%.2f, median=%.2f, quartiles=[%.2f-%.2f], range=[%d-%d]\n",mean(x),median(x),
                                              quantile(x,.25),quantile(x,.75),min(x),max(x)),sep="")
  x=as.vector(table(mdiffcal[,"cal"]))
  cat("   ---> num. scan per calib: ",sprintf("mean=%.2f, median=%.2f, quartiles=[%.2f-%.2f], range=[%d-%d]\n",mean(x),median(x),
                                              quantile(x,.25),quantile(x,.75),min(x),max(x)),sep="")
  
  invisible(list(match=mdiffcal,calib=unique(mdiffcal[,"calmz"]),calibori=lcalib))
}

.makeCalMatrices<-function(mdiffcal,lcalib,luscn){
  
  matmz=do.call("rbind",lapply(unique(mdiffcal[,"scan"]),function(ix){
    x=which(mdiffcal[,"scan"]==ix)
    tmp=tapply(mdiffcal[x,"mz",drop=F],mdiffcal[x,"calmz",drop=F],median)
    tmp[match(lcalib,names(tmp))]
  }))
  
  matint=do.call("rbind",lapply(unique(mdiffcal[,"scan"]),function(ix){
    x=which(mdiffcal[,"scan"]==ix)
    tmp=tapply(mdiffcal[x,"y",drop=F],mdiffcal[x,"calmz",drop=F],median)
    log10(tmp[match(lcalib,names(tmp))])
  }))
  matint=matint[match(luscn,unique(mdiffcal[,"scan"])),,drop=F]
  matmz=matmz[match(luscn,unique(mdiffcal[,"scan"])),,drop=F]
  mmppm=(sweep(matmz,2,lcalib,"/")-1)*10^6
  mmmz=(sweep(matmz,2,lcalib,"-"))
  dimnames(mmmz)=dimnames(mmppm)=dimnames(matint)=dimnames(matmz)=list(luscn,lcalib)
  
  return(list(mmmz,mmppm,matint,matmz))
  
}
chkCalib.plotSum<-function(mdiffcal,llspl=NULL,nsplit=1){
  
  lcalib=sort(unique(mdiffcal[,"calmz"]))
  if(is.null(llspl)) llspl=split(lcalib, ceiling(seq_along(1:length(lcalib))/nsplit))
  
  lcalib=unlist(llspl)
  
  luscn=as.character(min(mdiffcal[,1]):max(mdiffcal[,1]))

  re=.makeCalMatrices(mdiffcal,lcalib,luscn)
  mmmz=re[[1]]
  mmppm=re[[2]]
  matint=re[[3]]
  matmz=re[[4]]
  
  
  def.par <- par(no.readonly = TRUE) 
  par(mfrow=c(2,2),mar=c(3,6,1,.1),lwd=2,cex.main=1.2)
  
  for(type in 1:4){
 if(type==1){
   st1=apply(mmmz,2,quantile,c(.1,.25,.5,.75,.9),na.rm=T)
   xlab="Mass error [m/z]"
   xlabs=xlim=pretty(range(st1,na.rm=T))
   
 }
  
  if(type==2){
    st1=apply(mmppm,2,quantile,c(.1,.25,.5,.75,.9),na.rm=T)
    xlabs=xlim=pretty(range(st1,na.rm=T))
    xlab="Mass error [ppm]"
  }
  if(type==3){
    st1=apply(mmppm,2,function(x) quantile(which(!is.na(x)),c(.05,.25,.5,.75,.95),na.rm=T))
    xlabs=xlim=pretty(range(st1,na.rm=T))
    xlab="Measured scans"
  }
  
  if(type==4){
    st1=apply(matint,2,function(x) quantile(x,c(.1,.25,.5,.75,.9),na.rm=T))
    xlim=c(floor(min(st1,na.rm=T)),ceiling(max(st1,na.rm=T)))
    xlim=sort(c(log10(3*10^(xlim[1]:xlim[2]-1)),xlim[1]:xlim[2]))
    xlabs=sprintf("%0.e",10^xlim)
    
    xlab="Intensity"
  }
  
  imed=median(st1[3,])
  ylim=range(pretty(c(0,ncol(st1)*1.06)))
  plot(range(xlim),ylim,cex=0,xlab="",ylab="",axes=F,main=xlab)
  abline(h=rev(cumsum(sapply(llspl,length)))[-1]+.5,col="grey20",lty=3)
  segments(imed,0.5,imed,ncol(st1)+.5,lty=2,col="grey20",lwd=2)
  axis(1,at=xlim,labels = xlabs,pos=0)
  lines(st1[3,],1:ncol(st1),col="firebrick3")
  for(i in 1:ncol(st1)){
    segments(st1[1,i],i,st1[5,i],i,col="grey50")
    points(st1[,i],rep(i,5),pch=c(18,16,15,16,18),col=c("dodgerblue3","darkolivegreen3","firebrick3","darkolivegreen3","dodgerblue3"))
    axis(2,at=i,las=2,lab=colnames(st1)[i],tick=F,pos=xlim[1])
  }
  leg=c(rownames(st1)[2],paste(rownames(st1)[c(2,3)],collapse="/"),paste(rownames(st1)[c(1,5)],collapse="/"))
  legend("top",leg,pch=c(15,16,18),col=c("firebrick3","darkolivegreen3","dodgerblue3"),ncol=3,bty="n",pt.cex=1.2)
}
  par(def.par)
  
}


chkCalib.plotTrace<-function(mdiffcal,llspl=NULL,nsplit=1,what="dppm"){
  
  lcalib=sort(unique(mdiffcal[,"calmz"]))
  if(is.null(llspl)) llspl=split(lcalib, ceiling(seq_along(1:length(lcalib))/nsplit))

  
  def.par <- par(no.readonly = TRUE) 
  par(mar=c(3.4,3.4,1,.1),lwd=2,cex.main=1.2)
  
  ylim=pretty(mdiffcal[mdiffcal[,"calmz"]%in%unlist(llspl),what])
  
for(lcal in llspl){

tmpdat=mdiffcal[mdiffcal[,"calmz"]%in%lcal,]


xlim=pretty(range(tmpdat[,"scan"],na.rm=T))
colstit=brewer.pal(8,"Dark2")
plot(range(xlim),range(ylim),cex=0,xlab="",ylab="",
     ylim=range(c(0,ylim))+c(0,diff(range(ylim))/10),
     xlim=range(xlim)+c(0,diff(range(xlim))/5),
     axes=F,
     main=sprintf('Mass error [%s] for %.4f -> %.4f',what,min(lcal),max(lcal)))
segments(min(xlim),0,max(xlim),0,lwd=2)
for(i in seq(xlim[1],max(xlim),100)) segments(i,min(ylim),i,max(ylim),col="grey",lty=2,lwd=2)
axis(2,at=ylim,las=2,pos=xlim[1]);axis(1,at=xlim,pos=min(ylim))
#print(lcal)
imeds=do.call("rbind",tapply(tmpdat[,what],tmpdat[,"calmz"],quantile,c(.25,.5,.75),na.rm=T))
#print(imeds)
imeds=imeds[match(lcal,rownames(imeds)),,drop=F]
for(i in 1:nrow(imeds)){
  segments(min(xlim),imeds[i,1],max(xlim),imeds[i,1],col=colstit[i],lty=3,lwd=1)
  segments(min(xlim),imeds[i,3],max(xlim),imeds[i,3],col=colstit[i],lty=3,lwd=1)
  segments(min(xlim),imeds[i,2],max(xlim),imeds[i,2],col=colstit[i],lty=1,lwd=2)
}
for(i in 1:length(lcal)){
  l=which(tmpdat[,"calmz"]==lcal[i])
  lines(tmpdat[l,"scan"],tmpdat[l,what],lty=1,col=colstit[i])
  points(tmpdat[l,"scan"],tmpdat[l,what],pch=20,col=colstit[i])
}
leg=tapply(tmpdat[,1],tmpdat[,"calmz"],function(x) sprintf("N=%d [%d-%d]",length(x),min(x),max(x)))[as.character(lcal)]
legend(max(xlim),max(ylim),paste0(round(lcal,4),"\n",leg),pch=15,col=colstit,ncol=1,bty="n",pt.cex=1.2)
}
par(def.par)

}
