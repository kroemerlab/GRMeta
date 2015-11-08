.cvf <- function(x, n = 0) ifelse(sum(!is.na(x)) >= n, 100 * 
                                   sd(x, na.rm = T)/mean(x, na.rm = T), NA)
.cvrf <- function(x, n = 0) ifelse(sum(!is.na(x)) >= n, 100 * 
                                    mad(x, na.rm = T)/median(x, na.rm = T), NA)

corrSetQC<-function(obj,what,Samp2Corr=obj$Sid,Var2Corr=obj$Analyte,lQC=obj$Sid[which(obj$Meta$sType=="QC")],
                    nminQC=3,propNNA=0.5,lPcs=1:2,outfile=NULL,doplot=TRUE,Date2use="Date",complete="nothing",ipc=1,imod=4,verb=FALSE){
  
  if(!Date2use%in%names(obj$File)) stop("Date does not exist")
  Samp2Corr=Samp2Corr[!is.na(obj$File[Samp2Corr,Date2use])]
  m=obj$Data[[what]][Samp2Corr,Var2Corr]
  dts=obj$File[Samp2Corr,Date2use]
  dtsn=as.numeric(dts)#-mean(as.numeric(dts),na.rm=T)
  names(dtsn)=names(dts)=Samp2Corr
  
  
  clqc=lQC[lQC%in%Samp2Corr]
  clqc=clqc[order(dts[clqc])]
  curmat=names(which(colSums(!is.na(m[clqc,]))>=nminQC))
  clqc=clqc[(rowSums(!is.na(m[clqc,curmat]))/length(curmat))>propNNA]
  cat("QC samples:",clqc,"\n")
  cat("Num. of variables:",length(curmat),"out of",ncol(m),"\n")
  mqc=mqc2=log2(m[clqc,curmat])
  meqc=colMeans(mqc,na.rm=T)
  
  if(any(is.na(mqc2))) mqc2=t(impute.knn(t(mqc2),k=max(3,length(clqc)-2))$data)
  print(str(mqc))
  print(str(mqc2))
  mqc2=sweep(mqc2,2,meqc)
  pcqc=prcomp(mqc2)
  lPcs=lPcs[lPcs%in%(1:ncol(pcqc$rotation))]
  
  de=mean(dtsn[clqc])
  dts2=as.numeric(dts[clqc])-de
  andf=andfsa=list()
  for(i in lPcs){
    cat(".")
    ndf=data.frame(dts2=seq(min(dtsn),max(dtsn),length.out = 100)-de)
    ndfsa=data.frame(dts2=dtsn-de)
    
    ndfsa$dts2[ndfsa$dts2<min(ndf$dts2)]=min(ndf$dts2)
    ndfsa$dts2[ndfsa$dts2>max(ndf$dts2)]=max(ndf$dts2)
    
    lmm=lm(pcqc$x[,i]~dts2)
    ndf$pr1=predict(lmm,newdata=ndf,se=TRUE)$fit
    ndf$se1=predict(lmm,newdata=ndf,se=TRUE)$se.fit
    ndfsa$pr1=predict(lmm,newdata=ndfsa,se=TRUE)$fit
    ndfsa$se1=predict(lmm,newdata=ndfsa,se=TRUE)$se.fit
    ###
    lmm2=lm(pcqc$x[,i]~dts2+I(dts2^2))
    ndf$pr2=predict(lmm2,newdata=ndf,se=TRUE)$fit
    ndf$se2=predict(lmm2,newdata=ndf,se=TRUE)$se.fit
    ndfsa$pr2=predict(lmm2,newdata=ndfsa,se=TRUE)$fit
    ndfsa$se2=predict(lmm2,newdata=ndfsa,se=TRUE)$se.fit
    ###
    lmm3=lm(pcqc$x[,i]~dts2+I(dts2^2)+I(dts2^3))
    ndf$pr3=predict(lmm3,newdata=ndf,se=TRUE)$fit
    ndf$se3=predict(lmm3,newdata=ndf,se=TRUE)$se.fit
    ndfsa$pr3=predict(lmm3,newdata=ndfsa,se=TRUE)$fit
    ndfsa$se3=predict(lmm3,newdata=ndfsa,se=TRUE)$se.fit
    ###
    lmm4<-try(gam(pcqc$x[,i]~s(dts2)),TRUE)
    if("try-error"%in%class(lmm4)) lmm4=gam(pcqc$x[,i]~dts2)
    ndf$pr4=predict(lmm4,newdata=ndf,se=TRUE)$fit
    ndf$se4=predict(lmm4,newdata=ndf,se=TRUE)$se.fit
    ndfsa$pr4=predict(lmm4,newdata=ndfsa,se=TRUE)$fit
    ndfsa$se4=predict(lmm4,newdata=ndfsa,se=TRUE)$se.fit
    
    ndf$dts2=ndf$dts2+de
    ndfsa$dts2=ndfsa$dts2+de
    ndfsa$dts=dts
    
    andfsa[[as.character(i)]]=ndfsa
    andf[[as.character(i)]]=ndf
    
  }
  cat("\n")
  
  if(doplot | !is.null(outfile)) .incorrplot(pcqc,andf,andfsa,outfile=outfile,ipc=ipc)
  if(is.null(complete)) return(invisible(pcqc))
  if(!is.null(ipc) & !is.null(imod)){
    cat("Using model ",imod," on PC",ipc," to correct ",what,"\n",sep="")
    m0=obj$Data[[what]]
    ml=ml2=mlc=mlc2=log2(obj$Data[[what]])
    mlc[Samp2Corr,curmat]=sweep(ml[Samp2Corr,curmat],2,meqc[curmat])
    mlc2[Samp2Corr,curmat]=mlc[Samp2Corr,curmat]-
      andfsa[[as.character(ipc)]][Samp2Corr,paste("pr",imod,sep="")]%*%t(pcqc$r[curmat,ipc,drop=F])
    ml2[Samp2Corr,curmat]=sweep(mlc2[Samp2Corr,curmat],2,meqc[curmat],"+")
    ml2=2^ml2
    print(str(ml2))
    if(all(rownames(ml2)%in%Samp2Corr) & all(colnames(ml2)%in%Var2Corr)) cat("All samples and analytes from the original set were adjusted")
    if(any(!rownames(ml2)%in%Samp2Corr) | any(!colnames(ml2)%in%Var2Corr)){
      lsa=which(!rownames(ml2)%in%Samp2Corr)
   #   print(lsa)
      lv=which((!colnames(ml2)%in%Var2Corr))
    if(complete=="NA"){
      cat("Adding NAs for:\n")
      if(any((!rownames(ml2)%in%Samp2Corr))){ml2[lsa,]=NA;if(verb) cat(" *Sid:",rownames(ml2)[lsa],"\n")}
      if(any(!colnames(ml2)%in%Var2Corr)){ml2[lv,]=NA;if(verb) cat(" *Analytes:",colnames(ml2)[lv],"\n")}
    }
    if(complete=="remove"){
      cat("Removing:\n")
      if(any((!rownames(ml2)%in%Samp2Corr))){if(verb) cat(" *Sid:",rownames(ml2)[lsa],"\n");ml2=ml2[-lsa,,drop=F];}
      if(any(!colnames(ml2)%in%Var2Corr)){if(verb) cat(" *Analytes:",colnames(ml2)[lv],"\n");ml2=ml2[,-lv,drop=F]}
    }
    if(complete=="nothing"){
        cat("Former data used for:\n")
        if(any((!rownames(ml2)%in%Samp2Corr))){if(verb) cat(" *Sid:",rownames(ml2)[lsa],"\n")}
        if(any(!colnames(ml2)%in%Var2Corr)){if(verb) cat(" *Analytes:",colnames(ml2)[lv],"\n")}
    }
    }
    invisible(ml2)
  } else{invisible(pcqc)}
}

.incorrplot<-function(pcqc,andf,andfsa,outfile=NULL,ipc=1){
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  
  sumpc=summary(pcqc)
  curmat=rownames(pcqc$rotation)
  clqc=rownames(pcqc$x)
  lPCs=names(andf)
  
  if(!is.null(outfile)) pdf(outfile,width=10,height=7)
  
  for(i in lPCs){
    
    ndf=andf[[i]]
    ndfsa=andfsa[[i]]
    i=as.numeric(i)
    main=paste("PC",i,"(",round(sumpc$importance[2,i]*100,1),"%)",ifelse(i==ipc,"***",""))
    dts=ndfsa$dts;names(dts)=rownames(ndfsa)
    
    
    par(mar=c(5,4,1,.1),cex=1,cex.axis=1,cex.lab=1,lwd=2,cex.main=1)
    nf=layout(rbind(c(1,1),c(2,3)),widths = c(1,1),heights =  c(1,1), T)
    
    ##############
    ## Loadings
    ylim=range(pcqc$rotation[,i])
    ylim=ylim+c(-1,1)*diff(ylim)*.1
    out=outlierTest(lm(pcqc$rotation[,i]~1),n.max = 50,cutoff = 0.1)
    cols=rep("grey30",length(curmat))
    if(any(out$signif)){
      l=which(curmat%in%names(out$p))
      cols[l]="red"
    }
    
    xlim=c(0,max(length(curmat)+1,50))
    par(mar=c(1,4,1,.1))
    plot(1:length(curmat),pcqc$rotation[,i],col=cols,cex=1,axes=F,bty="n",ylab=main,xlab="",main=paste("Loadings on QC:",main),pch=16,ylim=ylim,xlim=xlim)
    abline(h=0)
    yaxt=axTicks(2)
    yaxt=c(yaxt,range(yaxt)+c(-1,1)*diff(yaxt[1:2]))
    axis(2,pretty(yaxt),las=2)
    if(any(out$signif)) 
      for(j in l) text(j,pcqc$rotation[j,i],gsub("@","\n@",curmat[j]),cex=.6,pos=c(1,3)[(sign(pcqc$rotation[j,i])==1)+1])
    
    
    ## 
    xlim=range(dts)
    xlim=chron(as.numeric(xlim)+c(-1,1)*as.numeric(diff(xlim)*.1))
    ylim=range(pcqc$x[,i])
    ylim[2]=ylim[2]+diff(ylim)*.25
    
    ############
    # Scores
    par(mar=c(5,4,1,.1))
    plot(dts[clqc],pcqc$x[,i],typ="b",axes=F,bty="n",ylab=main,main=paste("Fitting on QC:",main),xlab="Time",pch=16,xlim=xlim,ylim=ylim)
    yaxt=axTicks(2)
    yaxt=c(yaxt,range(yaxt)+c(-1,1)*diff(yaxt[1:2]))
    axis(2,pretty(yaxt),las=2)
    xaxt=axTicks(1)
    xaxt=sort(c(xaxt,range(xaxt)+c(-1,1)*diff(xaxt[1:2])))
    names(xaxt)=paste(hours(chron(xaxt)),":",minutes(chron(xaxt)),sep="")
    axis(1,at=xaxt,labels = names(xaxt))
    legend("top",c("Lin.","Quad.","Ord.3","Spline"),bty="n",col=1:4,ncol=4,pch=15,cex = .6)
    
    
    for(j in 1:4){
      lines(ndf$dts2,ndf[,paste("pr",j,sep="")],col=j,lwd=2)
      lines(ndf$dts2,ndf[,paste("pr",j,sep="")]+ndf[,paste("se",j,sep="")],col=j,lty=3)
      lines(ndf$dts2,ndf[,paste("pr",j,sep="")]-ndf[,paste("se",j,sep="")],col=j,lty=3)
    }
    text(dts[clqc],pcqc$x[,i],clqc,cex=.6,col=1,pos=3)
    
    ##############
    ## Prediction
    plot(xlim,ylim,cex=0,axes=F,bty="n",ylab=main,main=paste("Prediction:",main),xlab="Time",pch=16,xlim=xlim,ylim=ylim)
    axis(2,pretty(yaxt),las=2)
    axis(1,at=xaxt,labels = names(xaxt))
    legend("top",c("Lin.","Quad.","Ord.3","Spline"),bty="n",col=1:4,ncol=4,pch=15,cex=.6)
    
    for(j in 1:4)  points(ndfsa$dts2,ndfsa[,paste("pr",j,sep="")],col=j,cex=.5,pch=16)
    points(dts[clqc],pcqc$x[,i],pch=16,col=brewer.pal(8,"Accent")[6],cex=1.5)
    text(dts[clqc],pcqc$x[,i],clqc,cex=.6,col=brewer.pal(8,"Accent")[6],pos=3)
    par(mfrow=c(1,1))
  }
  if(!is.null(outfile)) dev.off()
  
  par(def.par)  #-
}
