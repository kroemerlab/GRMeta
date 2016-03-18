gatherOneEIC<-function(matEIC,sampfile,root="./",outfile=NULL,serial=TRUE){
  
  lSids=sampfile$Sid[file.exists(sampfile$Out)]
  
  d0=proc.time()[3]
  dfeic=eicinfos=sampinfos=list()
  if(serial) cat("Gathering ",length(unique(matEIC$Id))," EICs")
  for(ix in 1:length(lSids)){
    load(sampfile[ix,]$Out)
    if(!exists('mred')) next
    if(serial) cat(ifelse(ix%%10==0,"X","."))
    iattr=attributes(mred)[c("file","rts")]
    tabeic=attr(mred,"eic")
    if(any(tabeic$Id%in%matEIC$Id)){
      leics=which(tabeic$Id%in%matEIC$Id)
      mred=mred[mred[,"eic"]%in%leics,,drop=F]
      if(nrow(mred)==0) next
      lma=match(tabeic$Id,matEIC$Id) ##
      mred[,"eic"]=lma[mred[,"eic"]]
      tabeic=tabeic[leics,]
      mred=cbind(mred,samp=ix)
      eicinfos[[lSids[ix]]]=cbind(tabeic[order(match(tabeic$Id,matEIC$Id)),],samp=ix)
      sampinfos[[lSids[ix]]]=iattr
      dfeic[[lSids[ix]]]=mred
      rm(list=c("mred","iattr","tabeic"))
    }
  }
  
  if(serial) cat(" OK...\n")
  dfeic=do.call("rbind",dfeic)
  if(!serial) top=paste(length(unique(dfeic[,'eic'])),"/",nrow(matEIC)," in ",length(unique(dfeic[,'samp']))," samp. - ",sep="")
  
  if(!is.null(outfile)) if(!is.na(outfile)) save(file=outfile,dfeic,matEIC,eicinfos,sampinfos)
  
  if(!is.null(outfile)) if(is.na(outfile)){
    adfeic=dfeic
    for(inam in unique(matEIC$GrpEic)){
      
      if(serial) cat(inam," ")
      dfeic=data.frame(adfeic[adfeic[,"eic"]%in%which(matEIC$GrpEic==inam),,drop=FALSE])
      if(nrow(dfeic)==0) next
      dfeic$samp=lSids[dfeic$samp]
      dfeic$eic=matEIC$Id[dfeic$eic]
      #### expand missing scans
      grpss=paste(dfeic$samp,dfeic$eic,sep=";;")
      alladd=list()
      for(iss in unique(grpss)){
        l=which(grpss==iss)
        if(length(l)<2) next
        tmp=dfeic[l,,drop=F]
        lnewsc=min(tmp[,"scan"]):max(tmp[,"scan"])
        lnewsc=lnewsc[!lnewsc%in%tmp[,"scan"]]
        if(length(lnewsc)>0)
          alladd[[iss]]=data.frame(scan=lnewsc,mz=NA,y=0,rt=sampinfos[[tmp$samp[1]]]$rts[as.character(lnewsc)],
                                   eic=tmp$eic[1],samp=tmp$samp[1],stringsAsFactors=F)
        
      }
      if(length(alladd)>0){
        dfeic=rbind(dfeic,do.call("rbind",alladd))
        rownames(dfeic)=NULL
        dfeic=dfeic[order(dfeic$samp,dfeic$eic,dfeic$scan),]
        rm('alladd')
      }
      #### add mz/rtcor
      dfeic$mzcor=dfeic$mz
      dfeic$rtcor=dfeic$rt
      for(isamp in unique(dfeic$samp)){
        l=which(dfeic$samp==isamp)
        tmp=eicinfos[[isamp]]
        lma=match(dfeic$eic[l],tmp$Id)
        dfeic$mzcor[l]=dfeic$mz[l]*(1+tmp$shppm[lma]*10^-6)
        dfeic$rtcor[l]=dfeic$rt[l]+tmp$shrt[lma]
      }
      attr(dfeic,"eic")<-matEIC[matEIC$GrpEic==inam,]
      outfile=paste(root,"/",inam,"-isoSet.rda",sep="")
      save(file=outfile,dfeic,matEIC,eicinfos,sampinfos)
      
      
    }
    dfeic=adfeic
  }
  
  d1=proc.time()[3]
  
  if(!serial){
    top=paste(c(top,paste(round(d1-d0,1)," secs",sep="")),collapse="")
    return(list(top,d1-d0))
  }
  cat("Done in ",round(d1-d0,1)," secs\n",sep="")
  invisible(dfeic)
}

##########################

.GRgatherMultiEICCl<-function(ix,tabeic,sampfile,root){
  gatherOneEIC(tabeic[tabeic$GrpEic%in%ix,],sampfile,root=root,outfile=NA,serial=FALSE)
}

##########################

gatherMultiEICs<-function(matfile,tabeic,root="./",outfile=NA,ncl=1,nsplit=20,serial=TRUE){
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
  }

    if(!"GrpEic"%in%names(tabeic)) tabeic$GrpEic=tabeic$Id
  nmols=length(unique(tabeic$GrpEic))
  ngrps=max(1,round(nmols/nsplit))
  lst=seq(1,nmols,length.out=ngrps+1)
  llcode=suppressWarnings(split(unique(tabeic$GrpEic),1:ngrps))

  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  if(ncl==1){
    cat(" on 1 processor\n",sep="")
    lts=list()
    for(k in 1:length(llcode)){
        lts[[k]]=gatherOneEIC(tabeic[tabeic$GrpEic%in%llcode[[k]],],sampfile=matfile,root=root,outfile=outfile,serial=FALSE)
    }
    d1=proc.time()[3]
    cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(llcode),1)," secs per iso group\n",sep="")
    return(invisible(lts))
  }
  require("snowfall")
  
  cat(" on ",ncl," processors\n",sep="")
  sfInit(parallel=TRUE, cpus=ncl, type='SOCK',slaveOutfile=paste('GathEic',format(Sys.time(), "%y-%d-%b-%H:%M"),'.log',sep=""))
  sfLibrary(GRMeta)
  #sfExport( "gatherOneEIC", local=TRUE )
  re=sfClusterApplyLB(llcode,.GRgatherMultiEICCl,tabeic=tabeic,sampfile=matfile,root=root)
  sfStop()
  d1=proc.time()[3]
  cat(sapply(re,function(x) x[[1]]),sep="\n")
  lts=sapply(re,function(x) x[[2]])
  cat("Completed at ",date()," -> took ",round(d1-d0,1)," secs ",round((d1-d0)/length(llcode),1),"/",
      round((d1-d0)/length(llcode),1),"/",round(mean(lts),1)," secs per iso group\n",sep="")
  
  return(invisible(re))
  
}
