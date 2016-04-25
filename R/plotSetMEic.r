
plotSetMEic<-function(obj,WhichRT="rtcor",groupCol=NULL,colorCol=NULL,doPDF=T,addfile="./",endfile="-Eic",repDots="-",cexEL=0.6,...){
  
  # WhichRT="rtcor";groupCol=NULL;colorCol=NULL;addfile="./";endfile="-Eic";repDots="-";cexEL=0.6;dots=list()
  # width=10;height=14
  
  dots<-list(...)
  
  if(is.null(obj$Eic)) stop("Object does not contain any EIC infos.")
  
  lanalytes=NULL
  for(k in names(obj$Eic$File))  if(!is.null(dots[[k]])){
    l2add=obj$Eic$File$Analyte[obj$Eic$File[,k]%in%dots[[k]]]
    l2add=l2add[!l2add%in%lanalytes]
    if(length(l2add)>0) lanalytes=c(lanalytes,l2add)
  }
  if(is.null(lanalytes)){
    lanalytes=obj$Eic$File$Analyte
  }
  
  lanalytes=lanalytes[lanalytes%in%obj$Eic$File$Analyte]
  lfiles=lfiles2=unique(obj$Eic$File$EicFile[obj$Eic$File$Analyte%in%lanalytes])
  if(!is.null(obj$Eic$Path)) lfiles2=paste(obj$Eic$Path,lfiles2,sep="")
  lexist=which(file.exists(lfiles2))
  lfiles=lfiles[lexist]
  lanalytes=lanalytes[lanalytes%in%obj$Eic$File$Analyte[obj$Eic$File$EicFile%in%lfiles]]
  lfiles=lfiles[lfiles%in%obj$Eic$File$EicFile]
  cat(length(lanalytes),"analytes from",length(lfiles),"files\n")
  if(length(lanalytes)==0) stop("Nothing found!")
  
  width = ifelse(is.null(dots$width),5*par()$mfrow[2],dots$width)
  height = ifelse(is.null(dots$height),5*par()$mfrow[1],dots$height)
  dots=dots[names(dots)%in%names(par())]
  
  
  grp=obj$Meta$sType
  if(!is.null(groupCol)) if(!is.na(groupCol)) grp=obj$Meta[,groupCol]
  grp=factor(grp)
  llsids0=tapply(obj$Sid,grp,c)
  
  cols0=brewer.pal(9,"Set1")[-6][as.numeric(factor(obj$Meta$sType))]
  if(!is.null(colorCol)) if(!is.na(colorCol)) if(colorCol%in%names(obj$Meta)) cols0=obj$Meta[,colorCol]
  cols0[is.na(cols0)]="black"
  names(cols0)=obj$Sid
  
  par.def=par(no.readonly = TRUE)
  
  ifeic=lfiles[1]
  #print(ifeic)
  for(ifeic in lfiles){
    
    l2plots=obj$Eic$File$Analyte[obj$Eic$File$EicFile==ifeic]
    l2plots=l2plots[l2plots%in%lanalytes]
    if(!is.null(obj$Eic$Path)) ifeic=paste(obj$Eic$Path,ifeic,sep="")
    cat("Found",ifeic)
    dfeic<-eicpk<-NULL
    load(ifeic)
    cat(": ")
    
    
    dfeic=dfeic[which(dfeic$eic%in%obj$Eic$File[l2plots,]$Eic),]
    
    lse=paste(dfeic$samp,dfeic$eic)
    ## before/after peak zone
    dfeic$InPk2=dfeic$InPk3=FALSE
    dfeic$InPk2[unlist(tapply(1:nrow(dfeic),lse,function(x) x[x<=x[which(dfeic$InPk[x])[1]]]))]=TRUE
    dfeic$InPk3[unlist(tapply(1:nrow(dfeic),lse,function(x) x[x>=x[rev(which(dfeic$InPk[x]))[1]]]))]=TRUE
    
    lx=round(seq(1,length(l2plots),length.out = 5)[2:4])
    # if(sepx[2]%in%names(idf)) llsids=tapply(idf$Sid,idf[,sepx[2]],unique)
    #for(ix in 1:length(l2plots)){
    for(ix in 1:length(l2plots)){
      cat(ifelse(ix%in%lx,"X","."))
      analyte=l2plots[ix]
      ieic=obj$Eic$File[analyte,]$Eic
      
      imeth=ifelse(is.null(obj$Eic$File[analyte,]$Method),
                   obj$Annot[analyte,]$Method,obj$Eic$File[analyte,]$Method)
      if(!is.null(obj$Eic$Samp)){
        convsids=obj$Eic$Samp[,c("Sid",paste("SidEic",imeth,sep="."))]
      llsids=lapply(llsids0,function(x) convsids[match(x,convsids[,1]),2])
      cols=cols0
      names(cols)=convsids[match(names(cols),convsids[,1]),2]
      
      } else {
        llsids=llsids0
        cols=cols0
      }
      
      
      whichrt=WhichRT
      if(is.null(whichrt))  whichrt="rtcor"
      if(!whichrt%in%names(dfeic)) whichrt="rt"
      ceic=dfeic[dfeic$eic==ieic & dfeic$samp%in%unlist(llsids),]
      ceic$cols=cols[ceic$samp]
      rtr=range(ceic[,whichrt])
      rtr=rtr+c(-.5,1)*0.1*diff(rtr)
      rtr[1]=max(0,rtr[1])
      rtr=range(pretty(seq(rtr[1],rtr[2],length.out = 7)))
      
      Mint=obj$Eic$File[analyte,]$Bl
      
      outfile=paste(addfile,gsub("\\.",repDots,analyte),endfile,".pdf",sep="")
      if(doPDF) pdf(file=outfile,width=width,height=height)
      par(dots)
      .GRplotEIC(ceic,whichrt,llsids,Mint=Mint,rtr=rtr,cexEL=cexEL)
      if(doPDF) dev.off()
      
    }
    cat("\n")
    rm(list=c('dfeic','eicpk'))
    
  } ## end of for ifeic
  on.exit(par(par.def))
}




.GRplotEIC<-function(ceic,whichrt,llsids=list(All=unique(ceic$Sid)),Mint=NA,rtr=NULL,cexEL=0.4){
  
  
  alphadd<-function(hex.color.list,alpha=0.5) sprintf("%s%02X", hex.color.list, floor(alpha * 256))
  
  if(is.null(rtr)){
    rtr=range(ceic[,whichrt])
    rtr=rtr+c(-.5,1)*0.1*diff(rtr)
    rtr[1]=max(0,rtr[1])
    rtr=range(pretty(seq(rtr[1],rtr[2],length.out = 7)))
  }
  
  cols=tapply(ceic$cols,ceic$samp,unique)
  if(any(!ceic$InPk)) ceic$cols[!ceic$InPk]=alphadd(ceic$cols[!ceic$InPk],.4)
  
  for(namsamp in names(llsids)){
    #   cat(namsamp,"\n")
    lisamp=llsids[[namsamp]]
    lisamp=lisamp[lisamp%in%unique(ceic$samp)]# sample ids in ori file to considered in 
    #    lisamp[order(rtmat[lisamp,"HEap.1"],na.last = FALSE)]
    if(length(lisamp)==0){
      plot(0:1,0:1,axes=F,xlab="",ylab="",cex=0)
      text(.5,.5,namsamp)
      next
    }
    rint=range(ceic$y[ceic$samp%in%lisamp])*c(0,1.1)
    plot(rtr,rint,xlab=paste("Retention time (",whichrt,")",sep=""),
         ylab="Intensity (cps)",cex=0,bty="l",axes=FALSE,xlim=rtr,ylim=c(-rint[2]/50,rint[2]))
    abline(h=Mint,lty=2,lwd=par("lwd")*2)
    for(ik in lisamp){
      l=which(ceic$samp==ik & ceic$InPk2)
      if(length(l)>0) points(ceic[l,whichrt],ceic$y[l],col=ceic$cols[l],typ="l")
      l=which(ceic$samp==ik & ceic$InPk3)
      if(length(l)>0) points(ceic[l,whichrt],ceic$y[l],col=ceic$cols[l],typ="l")
    }
    leg=NULL
    for(ik in lisamp){
      l=which(ceic$samp==ik & ceic$InPk)
      if(length(l)>0){
        points(ceic[l,whichrt],ceic$y[l],col=ceic$cols[l],typ="l",lwd=par("lwd")*1.5)
        ir=range(ceic[l,whichrt])
        segments(ir,-rint[2]/50,ir,rint[2]/50,col=cols[ik],lwd=2)
        leg[ik]=paste(ik, " [",round(ir[1],2),'-',round(ir[2],2),"]",sep="")
      }
    }
    
    
    
    legend("topright",leg,pch=15,col=cols[lisamp],bty="n",cex=cexEL)
    legend("topleft",namsamp,bty="n",cex=par("cex.main"))
    xaxt=axTicks(1) #pretty(seq(rtr[1],rtr[2],length.out = 7))
    # xaxt=c(min(xaxt)-diff(xaxt)[1],xaxt,max(xaxt)+rev(diff(xaxt))[1])
    yaxt=pretty(c(0,rint[2]))
    yaxt=c(yaxt,max(yaxt)+rev(diff(yaxt))[1])
    axis(1,at=xaxt,las=1)
    axis(2,at=yaxt,las=2,pos=min(xaxt))
  } 
  
  
}

