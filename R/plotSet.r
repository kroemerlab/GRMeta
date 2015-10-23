

 ### Color Id

plot.metaboSet<-function(obj,outfile=NULL,
                         lgraphs=list(c("RT~InjOrder"),c("Area~InjOrder",log="y"),c("Area~1",log="y"),c("Height~Area",log="xy")),
                         mfrow=c(2,2),colorCol=NULL,deltaRT=0.05,linking="QC",orderBPlots=c("sType","InjOrder"),cexBP=0.5,cexBX=.8,cexEL=0.4,cexPT=1.2,...){
#  .plotOneAnalyte(obj,analyte = i,lgraphs=lgraphs,mfrow=mfrow,colorCol=colorCol,deltaRT=deltaRT,linking=linking,orderBPlots=orderBPlots,cexBP=cexBP,cexBX=cexBX,cexEL=cexEL,... )
  
  mgraphs=t(sapply(lgraphs,function(x) strsplit(x[1],"~")[[1]][1:2]))
  mgraphs[grep("^[0-9]$",mgraphs[,2]),2]=NA
  ltyps=c(names(obj$Meta),names(obj$Data))
  l2keep=which(mgraphs[,1]%in%ltyps & (mgraphs[,2]%in%ltyps | is.na(mgraphs[,2])))
  if("Eic" %in% names(obj)) l2keep=sort(c(l2keep,which(mgraphs[,1]=="Eic")))
#  print(mgraphs)
  lgraphs=lgraphs[l2keep]
  mgraphs=na.omit(unique(as.vector(mgraphs[l2keep,])))
  
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
    lanalytes=obj$Analyte
  }
  
  ###########
  ## Exclude some analytes that have no information on Data to be plotted
  l1=mgraphs[mgraphs%in%names(obj$Data)]
  m1=do.call("cbind",lapply(l1,function(x) colSums(!is.na(obj$Data[[x]][,lanalytes,drop=F]))))
  lanalytes=names(which(rowSums(m1>1)>0))
 if(any(rowSums(m1>1)==0)) cat("Excluded:",names(which(rowSums(m1>1)==0)),"\n")
  ###
  
  if(!is.null(outfile)){
      if(msg) cat("Analytes could not be matched: all are exported to ",outfile,"\n",sep="")
    pdf(file=outfile,width = ifelse(is.null(dots$width),8,dots$width),height = ifelse(is.null(dots$height),6,dots$height))
  }
  if(msg & is.null(outfile)){
    cat("Analytes could not be matched, using",lanalytes[1],"\n",sep="")
    lanalytes=lanalytes[1]
  }
  for(i in unique(lanalytes))
    .plotOneAnalyte(obj,analyte = i,lgraphs=lgraphs,mfrow=mfrow,colorCol=colorCol,deltaRT=deltaRT,linking=linking,orderBPlots=orderBPlots,
                    cexBP=cexBP,cexBX=cexBX,cexEL=cexEL,cexPT=cexPT,... )
  
  if(!is.null(outfile)) dev.off()
  
}


.infctlimlog=function(v){
  v=v[v>0]
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


################################
.plotOneAnalyte<-function(obj,analyte=obj$Analyte[1],
                         lgraphs=list(c("RT~InjOrder"),c("Area~InjOrder",log="y"),c("Area~1"),c("Height~Area",log="xy")),
                         mfrow=c(2,2),colorCol=NULL,deltaRT=0.05,linking="QC",orderBPlots=c("sType","InjOrder"),cexBP=0.5,cexBX=0.8,cexEL=0.4,cexPT=1.2,...){

# dots=list();analyte=obj$Analyte[1];lgraphs=list(c("RT~InjOrder"),c("Area~Height",log="yx"),c("Height~1"),c("Height~Sid",log="xy"));mfrow=c(2,2);deltaRT=0.05;linking=NULL;orderBPlots="sType";cexBP=0.5
# obj$Meta$Grp=factor(obj$Meta$sType,levels=c("Sa","QC"))
#print(lgraphs)
dots<-list(...)
dots=dots[names(dots)%in%names(par())]

lparams=unique(unlist(lapply(lgraphs,function(x) strsplit(x[1],"[~\\|]")[[1]])))
#print(lparams)
#print(lparams)
lparams1=unique(c(lparams[lparams%in%names(obj$Meta)],orderBPlots[orderBPlots%in%names(obj$Meta)]))
lparams2=lparams[lparams%in%names(obj$Data)]

idf=obj$Meta[,unique(c("Sid","sType",lparams1))]
for(i in lparams2)
  if(sum(!is.na(obj$Data[[i]][,analyte]))>2)
    eval(parse(text=paste("idf$",i,"=obj$Data[[i]][,analyte]",sep="")))

idf$color=brewer.pal(9,"Set1")[-6][as.numeric(factor(idf$sType))]
if(!is.null(colorCol)) if(!is.na(colorCol)) if(colorCol%in%names(obj$Meta)) idf$color=obj$Meta[,colorCol]

idf$color[is.na(idf$color)]="black"
rownames(idf)=idf$Sid
#print(str(idf))
if("InjOrder"%in%names(idf)) idf=idf[order(idf$InjOrder),]

par.def=par(no.readonly = TRUE)
#print(par.def)
###############
# Set RT/InjOrder/MZ
medRT0=injlim=NULL
.getrtlim<-function(idf,irt="RT",deltaRT){
  nround=ceiling(abs(log10(0.05)))
  medRT0=median(idf[,irt],na.rm=TRUE)
  medRT=round(median(idf[,irt],na.rm=TRUE),nround)
  rtlim=c(medRT0+c(-1,1)*deltaRT)
#  print(c(medRT0,medRT,deltaRT,rtlim))
  rtlim=c(floor(rtlim[1]*10^nround),ceiling(rtlim[2]*10^nround))/10^nround
  rtlim=seq(rtlim[1],rtlim[2],deltaRT/5)
  ly=range(idf$RT,medRT+deltaRT/2,medRT-deltaRT/2,na.rm=TRUE)
  rtlim=rtlim[rtlim>=min(ly) & rtlim<=max(ly)]
  return(list(rtlim=rtlim,medRT0=medRT0,medRT=medRT))
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
  x=lgraphs[[iplot]]
  if("log"%in%names(lgraphs[[iplot]])){
    logs=x["log"]
    x=x[!names(x)%in%"log"]
  }

  whatx=strsplit(x,"~")[[1]][2]
  whaty=strsplit(x,"~")[[1]][1]
  ########################
  #### y-axis
  if(!whaty%in%c(names(idf),"Eic")){plot.new();next}
  ylim=NULL
  if(whaty=="Eic") ylim=c(0,Inf)
  if(whaty!="Eic"){
    vy=idf[,whaty]
  if(whaty=="RT" | grepl("^RT\\.",whaty)){
    rtlims=.getrtlim(idf,whaty,deltaRT);rtlim=rtlims$rtlim;medRT0=rtlims$medRT0;medRT=rtlims$medRT;
    ylim=rtlim;vy[which(vy>=max(rtlim))]=max(rtlim);vy[which(vy<=min(rtlim))]=min(rtlim)
  }
  if(whaty=="InjOrder") ylim=injlim
  if(is.null(ylim)){
    if(length(grep("y",logs))>0){
      vy[which(vy<=0)]=min(vy,na.rm=TRUE)/2
      ylim=.infctlimlog(vy)
    } else ylim=pretty(seq(min(vy,na.rm=TRUE),max(vy,na.rm=TRUE),length.out = 7))
  }
  idf$Y=vy
  if(sum(!is.na(vy))<2){plot.new();next}
  }
#  print(ylim)
  ########################
  #### x-axis
  if(whaty!="Eic"){
  if(is.na(whatx) | !whatx%in%names(idf)){
    whatx=NULL
    logs=gsub("x","",logs)
  }
  xlim=vx=NULL
  if(!is.null(whatx) & whaty!="Eic"){
    vx=idf[,whatx]
    if(whatx=="RT" | grepl("^RT\\.",whatx)){
      rtlims=.getrtlim(idf,whatx,deltaRT);rtlim=rtlims$rtlim;medRT0=rtlims$medRT0;medRT=rtlims$medRT;
      xlim=rtlim;vx[which(vx>=max(rtlim))]=max(rtlim);vx[which(vx<=min(rtlim))]=min(rtlim)
    }
    if(whatx=="InjOrder") xlim=injlim
    if(is.character(vx)) xlim=c(0.3,length(vx)+.7)      
    if(is.factor(vx)) xlim=c(0.3,nlevels(vx)+.7)      
    if(is.null(xlim)){
      if(length(grep("x",logs))>0){
        vx[which(vx<=0)]=min(vx,na.rm=T)/2
        xlim=.infctlimlog(vx)
      } else xlim=pretty(seq(min(vx),max(vx),length.out = 7))
    }
    idf$X=vx
    if(!is.null(whatx)) if(sum(!is.na(idf[,whatx]))<2){plot.new();next}
  }
  }
  ########################
  par(dots)
#  print(dots)
  lsoSa=1:nrow(idf)
  sortSample=orderBPlots[orderBPlots%in%names(idf)]
  if(!is.null(sortSample)){
    lsoSa=paste("order(",paste(paste("idf$",sortSample,sep=""),collapse=","),")",sep="")
    lsoSa=eval(parse(text=lsoSa))
  }
  
  if(!is.null(whatx)){
    if(is.numeric(idf$X))
      .plotXY(idf,whatx,whaty,logs,xlim,ylim,analyte,linking,medRT0,deltaRT,cexPT)
    if(is.character(idf$X))
      .plotLinP(idf[lsoSa,],whaty,logs=gsub("x","",logs),xlim,ylim,analyte,cexBP,cexPT)
    if(is.factor(idf$X))
      .plotBoxP(idf[lsoSa,],whaty,logs=gsub("x","",logs),xlim,ylim,analyte,cexBX)

  }  
  if(is.null(whatx) & whaty!="Eic")
    .plotBarP(idf[lsoSa,],whaty,gsub("x","",logs),ylim,analyte,cexBP)

  ################
  if(whaty=="Eic"){
    
    sepx=strsplit(whatx,"\\|")[[1]]
    whichrt=sepx[1]
    llsids=list(All=unique(idf$Sid))
    if(sepx[2]%in%names(idf)) llsids=tapply(idf$Sid,idf[,sepx[2]],unique)

    ifeic=obj$Eic$File[analyte,]$EicFile
    if(is.null(obj$Eic$Path)) ifeic=paste(obj$Eic$Path,ifeic,sep="")
    ieicpk=obj$Eic$File[analyte,]$EicPK
    ieic=obj$Eic$File[analyte,]$EicId
    if(!file.exists(ifeic) | is.null(ieicpk) | is.na(ieicpk) | is.na(ieic)){
      for(i in names(llsids)){
        plot(0:1,0:1,axes=F,xlab="",ylab="",cex=0)
        text(.5,.5,i)
        next
      }
    }
    cat("Found",ifeic)
    load(ifeic)
    cat(".\n")
    
    if(is.null(whichrt))  whichrt="rtcor"
    if(whichrt%in%names(dfeic)) whichrt="rt"
    
    ceic=dfeic[dfeic$Eic==ieic & dfeic$Samp%in%obj$Eic$Sample$Samp,]
    ceic$Sid=obj$Eic$Sample$Sid[match(ceic$Samp,obj$Eic$Sample$Samp)]
    ceic$cols=idf[ceic$Sid,]$color
    
    rtr=range(ceic[,whichrt])
    rtr=rtr+c(-.5,1)*0.1*diff(rtr)
    rtr[1]=max(0,rtr[1])

    ipkmat=eicpk[[ieic]]$Pks[[ieicpk]]
    ipkmat=ipkmat[rownames(ipkmat)%in%obj$Eic$Sample$Samp,]
    rownames(ipkmat)=obj$Eic$Sample$Sid[match(rownames(ipkmat),obj$Eic$Sample$Samp)]
    Mint=obj$Eic$File[analyte,]$Mint
    .plotEIC(ceic,ipkmat,whichrt,llsids,Mint=Mint,rtr=rtr,cexEL=cexEL)
  } ### end of Eic
  
}

on.exit(par(par.def))
}

#######################################################################
#########################################
##### whatx is null
.plotBarP<-function(idf,whaty,logs="",ylim,analyte,cexBP){
  
  re=barplot(idf$Y,axes=F,ylab=whaty,bty="n",log=logs,ylim=range(ylim),main=analyte,col=idf$color,xpd=FALSE)
  axis(2,at=ylim,las=2)
  for(i in 1:nrow(re)) axis(1,at=re[i,1],labels = idf$Sid[i],cex.axis=cexBP,las=2,tick=F,pos=ylim[1])
}

.plotLinP<-function(idf,whaty,logs="",xlim,ylim,analyte,cexBP,cexPT=1.2){
  l=which(!is.na(idf$X))
  ylim2=which(ylim>min(idf$Y[l]) & ylim<max(idf$Y[l]))
  ylim2=max(min(ylim2)-1,1):min(max(ylim2)+1,length(ylim))
  
  re=plot(1:length(idf$X),idf$Y,axes=F,xlab="",ylab=whaty,bty="n",log=logs,ylim=range(ylim[ylim2]),xlim=range(xlim),main=analyte,col="grey20",type="l")
  points(1:length(idf$X),idf$Y,col=idf$color,pch=16,cex=cexPT)
  axis(2,at=ylim[ylim2],las=2)
  for(i in 1:nrow(idf)) axis(1,at=i,labels = idf$X[i],cex.axis=cexBP,las=2,tick=F,pos=min(ylim[ylim2]))
}

.plotBoxP<-function(idf,whaty,logs="",xlim,ylim,analyte,cexBX){
#  cat("Ylim: ",ylim)
  l=which(!is.na(idf$X))
  ylim2=which(ylim>min(idf$Y[l],na.rm=T) & ylim<max(idf$Y[l],na.rm=T))
  if(length(ylim2)<2) ylim2=1:length(ylim)
  ylim2=max(min(ylim2)-1,1):min(max(ylim2)+1,length(ylim))
  re=boxplot(Y~X,data=idf,axes=F,xlab="",ylab=whaty,bty="n",log=logs,ylim=range(ylim[ylim2]),xlim=range(xlim),main=analyte,cex=0)
  beeswarm(Y~X,data=idf,add=T,pwcol = idf$color,pch=16)
  axis(2,at=ylim[ylim2],las=2)
  labs=paste(re$names,"\n(",re$n,")",sep="")
  for(i in 1:length(labs)) axis(1,at=i,labels =labs[i],tick=F,cex.axis=cexBX)
}


##### whatx is numeric
.plotXY<-function(idf,whatx,whaty,logs="",xlim,ylim,analyte,linking=NULL,medRT0=NULL,deltaRT=NULL,cexPT=1.2){
  
  vx=idf$X
  vy=idf$Y
  
  plot(range(xlim),range(ylim),cex=0,axes=F,xlab=whatx,ylab=whaty,bty="n",log=logs,xlim=range(xlim),ylim=range(ylim),main=analyte)
  if(whatx=="RT" | grepl("^RT\\.",whatx)) abline(v=medRT0+c(-.5,0,.5)*deltaRT,lty=c(2,1,2),lwd=par("lwd")*2)
  if(whaty=="RT" | grepl("^RT\\.",whaty)) abline(h=medRT0+c(-.5,0,.5)*deltaRT,lty=c(2,1,2),lwd=par("lwd")*2)
  if(whatx=="InjOrder") abline(v=xlim,col="grey")
  if(whaty=="InjOrder") abline(h=ylim,col="grey")
  axis(1,at=xlim)
  axis(2,at=ylim,las=2)
  
  #if(!all(c(whatx,whaty)%in%c("RT","InjOrder"))){
  ndf=data.frame(x=seq(min(xlim),max(xlim),length.out = 100))
  trdf=data.frame(x=vx,y=vy)
  inlm<-function(form,form2,trdf,ndf){
    mod<-try(gam(form,data=trdf),TRUE)
    if( "try-error"%in%class(mod)) mod<-try(mod<-lm(form2,data=trdf),TRUE)
    if( "try-error"%in%class(mod)) return(NULL)
    if(any(as.character(form)=="log10(y)")) myfct<-function(i) 10^i else  myfct<-function(i) i
    pr=predict(mod,ndf,se.fit = T)
    pr=data.frame(pr)
    pr$x=ndf$x
    pr$fitu=myfct(pr$fit+1.96*pr$se)
    pr$fitl=myfct(pr$fit-1.96*pr$se)
    pr$fit=myfct(pr$fit)
    pr
  }
  
  if(logs=="") pr=inlm(y~s(x),y~x,trdf,ndf)
  if(logs%in%c("yx","xy"))  pr=inlm(log10(y)~s(log10(x)),log10(y)~log10(x),trdf,ndf)
  if(logs=="y")  pr=inlm(log10(y)~s(x),log10(y)~x,trdf,ndf)
  if(logs=="x")  pr=inlm(y~s(log10(x)),y~log10(x),trdf,ndf)
  
  if(!is.null(pr)){
  lines(pr$x,pr$fit,col="grey20")
  lines(pr$x,pr$fitl,col="grey30",lty=2)
  lines(pr$x,pr$fitu,col="grey30",lty=2)
  }
  points(vx,vy,col=idf$color,pch=16,cex=cexPT)
  if(!is.null(linking) & "InjOrder"==whatx)
    for(i in linking) lines(vx[which(idf$sType==i)],vy[which(idf$sType==i)],col=idf$color[which(idf$sType==i)][1])
  
}

############################
##
.plotEIC<-function(ceic,ipkmat,whichrt,llsids=list(All=unique(ceic$Sid)),Mint=NA,rtr=NULL,cexEL=0.4){
  
  
  alphadd<-function(hex.color.list,alpha=0.5) sprintf("%s%02X", hex.color.list, floor(alpha * 256))
  
  ceic$IsPk=ceic$IsPk2=ceic$IsPk3=FALSE
  for(i in rownames(ipkmat)) ceic$IsPk[which(ceic$Sid==i & ceic[,whichrt]>=ipkmat[i,"RTmi"] & ceic[,whichrt]<=ipkmat[i,"RTma"])]=TRUE
  for(i in rownames(ipkmat)) ceic$IsPk2[which(ceic$Sid==i & (ceic[,whichrt]<=ipkmat[i,"RTmi"]))]=TRUE
  for(i in rownames(ipkmat)) ceic$IsPk3[which(ceic$Sid==i & (ceic[,whichrt]>=ipkmat[i,"RTma"]))]=TRUE
  
  if(is.null(rtr)){
    rtr=range(ceic[,whichrt])
    rtr=rtr+c(-.5,1)*0.1*diff(rtr)
    rtr[1]=max(0,rtr[1])
    rtr=range(pretty(seq(rtr[1],rtr[2],length.out = 7)))
  }
  
  rtmat=ipkmat[,c("Samp","RTap","RTap.1","RTmi","RTma","HEap.1","HEap")]
  cols=tapply(ceic$cols,ceic$Sid,unique)
  rtmat$cols=cols[rownames(ipkmat)]
  
  for(namsamp in names(llsids)){
    #   cat(namsamp,"\n")
    lisamp=llsids[[namsamp]]
    lisamp=lisamp[lisamp%in%unique(ceic$Sid)]# sample ids in ori file to considered in 
    lisamp[order(rtmat[lisamp,"HEap.1"],na.last = FALSE)]
    if(length(lisamp)==0){
      plot(0:1,0:1,axes=F,xlab="",ylab="",cex=0)
      text(.5,.5,namsamp)
      next
    }
    rint=range(ceic$y[ceic$Sid%in%lisamp])*c(0,1.1)
    plot(rtr,rint,xlab=paste("Retention time (",whichrt,")",sep=""),ylab="Intensity (cps)",cex=0,bty="l",axes=FALSE,xlim=rtr)
    abline(h=Mint,lty=2,lwd=par("lwd")*2)
    for(ik in lisamp){
      points(ceic[,whichrt][ceic$Sid==ik & ceic$IsPk2],ceic$y[ceic$Sid==ik & ceic$IsPk2],col=alphadd(cols[ik],.4),typ="l")
      points(ceic[,whichrt][ceic$Sid==ik & ceic$IsPk3],ceic$y[ceic$Sid==ik & ceic$IsPk3],col=alphadd(cols[ik],.4),typ="l")
    }
    for(ik in lisamp)
      points(ceic[,whichrt][ceic$Sid==ik & ceic$IsPk],ceic$y[ceic$Sid==ik & ceic$IsPk],col=cols[ik],typ="l",lwd=par("lwd")*1.5)
    
    toprt=apply(rtmat[lisamp,c("RTmi","RTma")],1,function(x) sprintf("%.2f-%.2f",x[1],x[2]))
    toprt[which(toprt=="NA-NA")]=""
    legend("topright",paste(lisamp,toprt),pch=15,col=cols[lisamp],bty="n",cex=cexEL)
    lma3=match(lisamp,rtmat$Samp)
    points(rtmat[lisamp,]$RTap.1,rtmat[lisamp,]$HEap.1,col=rtmat[lisamp,]$cols,pch=17)
    #      points(rtmat$RTap[lma3],rtmat$HEap[lma3],col=meta$cols[lma2],pch=16,cex=paramsP$cexpt)
    segh=rint[2]/50
    segments(rtmat[lisamp,]$RTma,-segh,rtmat[lisamp,]$RTma,segh,col=rtmat[lisamp,]$cols,lwd=par("lwd")*2)
    segments(rtmat[lisamp,]$RTmi,-segh,rtmat[lisamp,]$RTmi,segh,col=rtmat[lisamp,]$cols,lwd=par("lwd")*2)
    
    legend("topleft",namsamp,bty="n",cex=par("cex.main"))
    xaxt=axTicks(1) #pretty(seq(rtr[1],rtr[2],length.out = 7))
   # xaxt=c(min(xaxt)-diff(xaxt)[1],xaxt,max(xaxt)+rev(diff(xaxt))[1])
    yaxt=axTicks(2)
    yaxt=c(yaxt,max(yaxt)+rev(diff(yaxt))[1])
    axis(1,at=xaxt,las=1)
    axis(2,at=yaxt,las=2,pos=min(xaxt))
  } 
  
  
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
