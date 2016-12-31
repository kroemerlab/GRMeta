

.inChkCalibPlot1<-function(mdiffcal,llspl){
  
  lcalib=unlist(llspl)
  luscn=as.character(min(mdiffcal[,1]):max(mdiffcal[,1]))
  matmz=do.call("rbind",tapply(1:nrow(mdiffcal),mdiffcal[,"scan"],function(x){
    tmp=mdiffcal[x,]
    tmp[match(lcalib,tmp[,"calmz"]),"mz"]
  }))
  matint=do.call("rbind",tapply(1:nrow(mdiffcal),mdiffcal[,"scan"],function(x){
    tmp=mdiffcal[x,]
    log10(tmp[match(lcalib,tmp[,"calmz"]),"y"])
  }))
  matint=matint[match(luscn,rownames(matint)),,drop=F]
  matmz=matmz[match(luscn,rownames(matmz)),,drop=F]
  mmppm=(sweep(matmz,2,lcalib,"/")-1)*10^6
  mmmz=(sweep(matmz,2,lcalib,"-"))
  dimnames(mmmz)=dimnames(mmppm)=dimnames(matint)=dimnames(matmz)=list(luscn,lcalib)
  
  
  
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


.inChkCalibPlot2<-function(mdiffcal,llspl,what="dppm"){
  
  
  def.par <- par(no.readonly = TRUE) 
  par(mar=c(3.4,3.4,1,.1),lwd=2,cex.main=1.2)
  
  ylim=pretty(mdiffcal[mdiffcal[,"calmz"]%in%unlist(llspl),what])
  
for(lcal in llspl){

tmpdat=mdiffcal[mdiffcal[,"calmz"]%in%lcal,]
what="dppm"

xlim=pretty(range(tmpdat[,"scan"],na.rm=T))
colstit=brewer.pal(8,"Dark2")
plot(range(xlim),range(ylim),cex=0,xlab="",ylab="",ylim=range(ylim)+c(0,diff(range(ylim))/10),axes=F,
     main=sprintf('Mass error [%s] for %.4f -> %.4f',what,min(lcal),max(lcal)))
abline(h=0)
for(i in seq(xlim[1],max(xlim),100)) segments(i,min(ylim),i,max(ylim),col="grey",lty=2,lwd=2)
axis(2,at=ylim,las=2,pos=xlim[1]);axis(1,pos=min(ylim))
imeds=do.call("rbind",tapply(tmpdat[,what],tmpdat[,"calmz"],quantile,c(.25,.5,.75),na.rm=T))
imeds=imeds[match(lcal,rownames(imeds)),]
abline(h=imeds[,1],col=colstit,lty=3,lwd=1)
abline(h=imeds[,3],col=colstit,lty=3,lwd=1)
abline(h=imeds[,2],col=colstit,lty=2)
for(i in 1:length(lcal)){
  l=which(tmpdat[,"calmz"]==lcal[i])
  lines(tmpdat[l,"scan"],tmpdat[l,what],lty=1,col=colstit[i])
  points(tmpdat[l,"scan"],tmpdat[l,what],pch=20,col=colstit[i])
}
leg=tapply(tmpdat[,1],tmpdat[,"calmz"],function(x) sprintf(" (N=%d, sc.=%d-%d)",length(x),min(x),max(x)))[as.character(lcal)]
legend("top",paste0(round(lcal,4),leg),pch=15,col=colstit,ncol=3,bty="n",pt.cex=1.2)
}
par(def.par)

}
