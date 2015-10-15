plot.metaboSet<-function(obj,outfile=NULL,lgraphs=list(c("InjOrder","RT"),c("InjOrder","Area"),c("Area"),c("Height","Area")),
                mfrow=c(2,2),deltaRT=0.05,linking="QC",orderBPlots="sType",cexBP=0.5,...){
  
  lgraphs=lgraphs[sapply(lgraphs,function(x) all(x%in%c(names(obj$Meta),names(obj$Data))))]
  if(length(lgraphs)==0) stop("Graphs cannot be matched to the obj$Meta or obj$Data")

  dots<-list(...)
  lanalytes=NULL
  for(k in names(obj$Annot))  if(!is.null(dots[[k]])){
    l2add=obj$Annot$Analyte[obj$Annot[,k]%in%dots[[k]]]
    l2add=l2add[!l2add%in%lanalytes]
    if(length(l2add)>0) lanalytes=c(lanalytes,l2add)
  }
  
  msg=FALSE
  if(is.null(lanalytes)){
    msg=TRUE
    lanalytes=obj$Annot$Analyte
  }
  
  ###########
  l1=unique(unlist(lgraphs)[unlist(lgraphs)%in%names(obj$Data)])
  m1=do.call("cbind",lapply(l1,function(x) colSums(!is.na(obj$Data[[x]][,lanalytes,drop=F]))))
  lanalytes=names(which(rowSums(m1>2)>0))

  if(!is.null(outfile)){
      if(msg) cat("Analytes could not be matched: all are exported to ",outfile,"\n",sep="")
    pdf(file=outfile,width = ifelse(is.null(dots$width),8,dots$width),height = ifelse(is.null(dots$height),6,dots$height))
  }
  if(msg & is.null(outfile)){
    cat("Analytes could not be matched, using",lanalytes[1],"\n",sep="")
    lanalytes=lanalytes[1]
  }
  for(i in unique(lanalytes))
    .plotOneAnalyte(obj,analyte = i,lgraphs=lgraphs,mfrow=mfrow,deltaRT=deltaRT,linking=linking,orderBPlots=orderBPlots,... )
  
  if(!is.null(outfile)) dev.off()
  
}





.plotOneAnalyte<-function(obj,analyte=obj$Analyte[1],
                         lgraphs=list(c("InjOrder","RT"),c("InjOrder","Area"),c("Height"),c("Height","Area")),
                         mfrow=c(2,2),deltaRT=0.05,linking="QC",orderBPlots="sType",...){

infctlim=function(v){
  lly=rep(c(1,2,5),31)*10^rep(-10:20,each=3)
  llyt=as.character(lly)
  names(lly)=llyt
  
  lly2=rep(c(1),31)*10^rep(-10:20,each=1)
  lly2t=as.character(lly2)
  names(lly2)=lly2t
  
  helim=lly[sort(min(which(lly>min(v,na.rm=T))):max(which(lly<max(v,na.rm=T))))]
  helim=c(max(lly[lly<helim[1]]),helim,min(lly[lly>max(helim)]))
  if(length(helim)>7){
    helim=lly2[sort(min(which(lly2>min(v,na.rm=T))):max(which(lly2<max(v,na.rm=T))))]
    helim=c(max(lly2[lly2<helim[1]]),helim,min(lly2[lly2>max(helim)]))
  }
  return(helim)
}


dots<-list(...)
dots=dots[names(dots)%in%names(par())]

lparams=unique(unlist(lgraphs))
lparams1=lparams[lparams%in%names(obj$Meta)]
lparams2=lparams[lparams%in%names(obj$Data)]

idf=obj$Meta[,c("Sid","sType",lparams1)]
for(i in lparams2)
  if(sum(!is.na(obj$Data[[i]][,analyte]))>2)
    eval(parse(text=paste("idf$",i,"=obj$Data[[i]][,analyte]",sep="")))

if(!is.null(obj$Meta$Color)){idf$color=obj$Meta$Color}else{
  idf$color=brewer.pal(9,"Set1")[-6][as.numeric(factor(idf$sType))]
  idf$color[is.na(idf$color)]="black"
}
idf=idf[order(idf$InjOrder),]

par.def=par(no.readonly = TRUE)
###############
# source('~/Metabo/Aspirin3/qqq/plotfct.r')

if("RT"%in%lparams){
  nround=ceiling(abs(log10(0.05)))
  medRT0=median(idf$RT,na.rm=TRUE)
  medRT=round(median(idf$RT,na.rm=TRUE),nround)
  rtlim=c(medRT0+c(-1,1)*deltaRT)
  rtlim=c(floor(rtlim[1]*10^nround),ceiling(rtlim[2]*10^nround))/10^nround
  rtlim=seq(rtlim[1],rtlim[2],deltaRT/5)
  ly=range(idf$RT,medRT+deltaRT/2,medRT-deltaRT/2,na.rm=TRUE)
  rtlim=rtlim[rtlim>=min(ly) & rtlim<=max(ly)]
}

if("InjOrder"%in%lparams){
  injlim=range(idf$InjOrder,na.rm=TRUE)*c(.9,1)
  f=ifelse(diff(injlim)>80,20,10)
  injlim=seq(floor(injlim[1]/f),ceiling(injlim[2]/f)*f,f)
}
#####

if(!is.null(mfrow)) par(mfrow=mfrow)
for(iplot in 1:length(lgraphs)){
  logs=""
  whatx=whaty=lgraphs[[iplot]][1]
  if(length(lgraphs[[iplot]])==2) whaty=lgraphs[[iplot]][2]
  
  if(!all(lgraphs[[iplot]]%in%names(idf))){plot.new();next}
  
  if(sum(!is.na(idf[,whatx]) & !is.na(idf[,whaty]))<2){plot.new();next}

  vx=idf[,whatx]
  vy=idf[,whaty]
  
if(whatx=="RT"){xlim=rtlim;vx[which(vx>=max(rtlim))]=max(rtlim);vx[which(vx<=min(rtlim))]=min(rtlim)}
if(whatx=="InjOrder") xlim=injlim
if(!whatx%in%c("InjOrder","RT")){vx[which(vx<=0)]=min(vx,na.rm=T);xlim=infctlim(vx);logs="x"}

if(whaty=="RT"){ylim=rtlim;vy[which(vy>=max(rtlim))]=max(rtlim);vy[which(vy<=min(rtlim))]=min(rtlim)}
if(whaty=="InjOrder") ylim=injlim
if(!whaty%in%c("InjOrder","RT")){vy[which(vy<=0)]=min(vy,na.rm=T);ylim=infctlim(vy);logs=paste(logs,"y",sep="")}

if(whatx!=whaty){
  par(dots)
  plot(range(xlim),range(ylim),cex=0,axes=F,xlab=whatx,ylab=whaty,bty="n",log=logs,xlim=range(xlim),ylim=range(ylim),main=analyte)
if(whatx=="RT") abline(v=medRT0+c(-.5,0,.5)*deltaRT,lty=c(2,1,2),lwd=par()$lwd*2)
if(whaty=="RT") abline(h=medRT0+c(-.5,0,.5)*deltaRT,lty=c(2,1,2),lwd=par()$lwd*2)
if(whatx=="InjOrder") abline(v=injlim,col="grey")
if(whaty=="InjOrder") abline(h=injlim,col="grey")
axis(1,at=xlim)
axis(2,at=ylim,las=2)
if(logs=="xy"){
  pr=predict(lm(log10(y)~log10(x),data=data.frame(x=vx,y=vy)),newdata = data.frame(x=xlim))
  lines(xlim,10^pr)
}
points(vx,vy,col=idf$color,pch=16)
if(!is.null(linking) & "InjOrder" %in%lgraphs[[iplot]])
  for(i in linking) lines(vx[which(idf$sType==i)],vy[which(idf$sType==i)],col=idf$color[which(idf$sType==i)][1])

} ## end of plot
  
  
  if(whatx==whaty){
    par(dots)
    lso=idf$Sid
    if(orderBPlots%in%names(obj$Meta)) lso=obj$Sid[order(obj$Meta[,orderBPlots],obj$Meta$InjOrder)]
    lso=match(lso,idf$Sid)
    re=barplot(vx[lso],axes=F,ylab=whatx,bty="n",log=ifelse(logs=="","","y"),ylim=range(xlim),main=analyte,col=idf$color[lso])
    axis(2,at=xlim,las=2)
    for(i in 1:nrow(re)) axis(1,at=re[i,1],labels = idf$Sid[lso][i],cex.axis=cexBP,las=2,tick=F,pos=xlim[1])
  }


}

par(par.def)

}

##################
# plotManyAnalyte<-function(obj,outfile=NULL,width=8,height=6, Analyte=obj$Analyte,MetName=NULL,Method=NULL,...){
#   
#   
#   lanalytes=NULL
#   if(!is.null(Analyte)) lanalytes=Analyte[Analyte%in%obj$Annot$Analyte]
#   if(!is.null(MetName)) lanalytes=c(lanalytes,obj$Annot$MetName[obj$Annot$MetName%in%MetName])
#   if(!is.null(Method)) lanalytes=c(lanalytes,obj$Annot$Method[obj$Annot$Method%in%Method])
#   if(!is.null(outfile)) pdf(file=outfile,width = width,height = 8)
#   for(i in lanalytes) plotOneAnalyte(obj,analyte = i,... )
#   if(!is.null(outfile)) dev.off()
#   
# }
# 
