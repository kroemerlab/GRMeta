gatherOneEIC<-function(matEIC,sampfile,outfile=NULL,eicParams,doMerge=TRUE,verbose=TRUE,serial=TRUE){
  
  lSids=sampfile$Sid[file.exists(sampfile$Out)]
  
  d0=proc.time()[3]
  dfeic=eicinfos=sampinfos=list()
  if(verbose) cat("Gathering ",length(unique(matEIC$Id))," EICs")
  for(ix in 1:length(lSids)){
    load(sampfile[ix,]$Out)
    if(!exists('mred')) next
    if(verbose) cat(ifelse(ix%%10==0,"X","."))
    iattr=attributes(mred)[c("file","rts","addRT","addMZ" )]
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
  
  if(verbose) cat(" OK...\n")
  dfeic=do.call("rbind",dfeic)
  if(!serial) top=paste(length(unique(dfeic[,'eic'])),"/",nrow(matEIC)," in ",length(unique(dfeic[,'samp']))," samp. - ",sep="")
  
  if(!is.null(outfile)) if(!is.na(outfile)) save(file=outfile,dfeic,matEIC,eicinfos,sampinfos)
  
  if(!is.null(outfile)) if(is.na(outfile)){
    adfeic=dfeic
    for(inam in unique(matEIC$GrpEic)){
      
      if(verbose) cat(inam," ")
      dfeic=data.frame(adfeic[adfeic[,"eic"]%in%which(matEIC$GrpEic==inam),,drop=FALSE])
      if(nrow(dfeic)==0) next
      dfeic$samp=lSids[dfeic$samp]
      dfeic$eic=matEIC$Id[dfeic$eic]
      # #### add mz/rtcor
      # dfeic$mzcor=dfeic$mz
      # dfeic$rtcor=dfeic$rt
      # for(isamp in unique(dfeic$samp)){
      #   l=which(dfeic$samp==isamp)
      #   tmp=eicinfos[[isamp]]
      #   lma=match(dfeic$eic[l],tmp$Id)
      #   dfeic$mzcor[l]=dfeic$mz[l]*(1+tmp$shppm[lma]*10^-6)
      #   dfeic$rtcor[l]=dfeic$rt[l]+tmp$shrt[lma]
      # }
      #### remerging 
      whichmz=ifelse("mzcor"%in%names(dfeic),"mzcor","mz")
      whichrt=ifelse("rtcor"%in%names(dfeic),"rtcor","rt")
      if(doMerge) dfeic=.GRmergingEIC(dfeic,eicParams,whichmz=whichmz, whichrt=whichrt,verbose=verbose)
      #### expand missing scans
      dfeic=.GRaddmiss(dfeic,sampinfos)
      dfeic=dfeic[order(dfeic$eic,dfeic$samp,dfeic$scan,dfeic$mzcor),]
      #### Compute eic stats
      l=which(dfeic$y>0)
      if("ineic" %in% names(dfeic)) l=l[which(dfeic[,"ineic"]==1)]
      eicst=data.frame(GrpEic=inam,Id=unname(tapply(dfeic$eic[l],dfeic$eic[l],unique)),
                 do.call("rbind",tapply(l,dfeic$eic[l],function(x) 
                   c(dfeic[x[which.max(dfeic$y[x])],whichrt],range(dfeic[x,whichrt])))),
                 do.call("rbind",tapply(l,dfeic$eic[l],function(x) 
                   c(dfeic[x[which.max(dfeic$y[x])],whichmz],range(dfeic[x,whichmz])))),
                 stringsAsFactors=F,row.names=NULL)
      names(eicst)=c("GrpEic", "Id","rtap","rtmin","rtmax","mzap","mzmin","mzmax" )
      attr(dfeic,"oeic")<-matEIC[matEIC$GrpEic==inam,]
      attr(dfeic,"eic")<-eicst
      outfile=paste(ifelse(is.null(eicParams$dirEic),".",eicParams$dirEic),"/",inam,eicParams$addEic,".rda",sep="")
      save(file=outfile,dfeic,eicinfos,sampinfos)
      
      #### end of for in each group
    }
    if(verbose) cat(inam,"\n")
    dfeic=adfeic
  }
  
  d1=proc.time()[3]
  
  if(!serial){
    top=paste(c(top,paste(round(d1-d0,1)," secs",sep="")),collapse="")
    return(list(top,d1-d0))
  }
  if(verbose) cat("Done in ",round(d1-d0,1)," secs\n",sep="")
  invisible(dfeic)
}

##########################

.GRgatherMultiEICCl<-function(ix,tabeic,sampfile,eicParams,doMerge){
  gatherOneEIC(tabeic[tabeic$GrpEic%in%ix,],sampfile,outfile=NA,eicParams=eicParams,doMerge=doMerge,verbose=FALSE,serial=FALSE)
}

##########################

gatherMultiEICs<-function(matfile,tabeic,outfile=NA,eicParams,doMerge=TRUE,ncl=1,nsplit=20,serial=TRUE,verbose=TRUE){
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,nrow(matfile),parallel:::detectCores()))
  }

    if(!"GrpEic"%in%names(tabeic)) tabeic$GrpEic=tabeic$Id

      llcode=split(unique(tabeic$GrpEic), ceiling(seq_along(unique(tabeic$GrpEic))/nsplit))
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  if(ncl==1){
    cat(" on 1 processor\n",sep="")
    lts=list()
    for(k in 1:length(llcode)){
        lts[[k]]=gatherOneEIC(tabeic[tabeic$GrpEic%in%llcode[[k]],],sampfile=matfile,outfile=outfile,eicParams=eicParams,doMerge=doMerge,serial=FALSE,verbose=verbose)
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
  re=sfClusterApplyLB(llcode,.GRgatherMultiEICCl,tabeic=tabeic,sampfile=matfile,eicParams=eicParams,doMerge=doMerge)
  sfStop()
  d1=proc.time()[3]
  cat(sapply(re,function(x) x[[1]]),sep="\n")
  lts=sapply(re,function(x) x[[2]])
  cat("Completed at ",date()," -> took ",round(d1-d0,1)," secs ",round((d1-d0)/length(llcode),1),"/",
      round((d1-d0)/length(llcode),1),"/",round(mean(lts),1)," secs per iso group\n",sep="")
  
  return(invisible(re))
  
}


## need sampinfos
.GRaddmiss<-function(df,sampinfos){
  grpss=paste(df$samp,df$eic,sep=";;")
  alladd=list()
  for(iss in unique(grpss)){
    l=which(grpss==iss)
    if(length(l)<2) next
    tmp=df[l,,drop=F]
    l2use=1:nrow(tmp)
    if("ineic" %in% names(df)) l2use=which(tmp[,"ineic"]==1)
    if(length(l2use)<2) next
    lnewsc=min(tmp[l2use,"scan"]):max(tmp[l2use,"scan"])
    lnewsc=lnewsc[!lnewsc%in%tmp[l2use,"scan"]]
    if(length(lnewsc)>0){
      alladd[[iss]]=data.frame(scan=lnewsc,mz=NA,y=0,rt=sampinfos[[tmp$samp[1]]]$rts[as.character(lnewsc)],
                               eic=tmp$eic[1],samp=tmp$samp[1],stringsAsFactors=F)
      if("rtcor"%in%names(tmp))  alladd[[iss]]$rtcor=alladd[[iss]]$rt+mean(tmp$rtcor-tmp$rt,na.rm=T)
      if("mzcor"%in%names(tmp))  alladd[[iss]]$mzcor=alladd[[iss]]$mz*mean(tmp$mzcor/tmp$mz,na.rm=T)
      if("ineic"%in%names(tmp))  alladd[[iss]]$ineic=1
    }
    
  }
  if(length(alladd)>0){
    alladd=do.call("rbind",alladd)
    df=rbind(df,alladd[,names(df)])
    rownames(df)=NULL
    df=df[order(df$samp,df$eic,df$scan),]
    rm('alladd')
  }
  return(df)
}

###
.GRmergeInterval<-function(m){
  if(is.null(rownames(m))) rownames(m)=1:nrow(m)
  le=rep(rownames(m),2)[order(as.vector(m))]
  len=rep(NA,length(le))
  len[which(le==le[1])]=1
  while(any(is.na(len))){
    j=which(is.na(len))[1]
    k=max(which(!is.na(len)))
    if(k>j) len[which(le==le[j])]=len[k]
    if(k<j) len[which(le==le[j])]=max(len,na.rm=T)+1
  }
  tapply(le,len,unique)
}

#########
.GRmergingEIC<-function(dfeic,eicParams,whichmz="mzcor",whichrt="rtcor",verbose=TRUE){
  
  if("ineic" %in% names(dfeic)) l=which(!is.na(dfeic[,whichmz]) & dfeic[,"ineic"]==1)
  if(!"ineic" %in% names(dfeic)) l=which(!is.na(dfeic[,whichmz]))

  lpp=.GRsplist(dfeic[l,whichmz],l,ismass=TRUE,d=eicParams$dPPM)
  ldups=lapply(lpp,function(x) unique(dfeic$eic[x]))
  l2excl=names(which(table(unlist(ldups))>1))
  if(length(l2excl)>0){ldups=lapply(ldups,function(x) x[!x%in%l2excl])
  ldups=c(ldups,as.list(l2excl))
  }
  ldups2=ldups
  ndups=sapply(ldups,length)
  ### 
  if(any(ndups>1)){
    for(i in which(ndups>1)){
  li=l[dfeic$eic[l]%in%ldups[[i]]]
  lirt=do.call("cbind",tapply(dfeic[li,whichrt],dfeic[li,"eic"],range))
  lirt=lirt[,order(lirt[1,]),drop=F]
  lcuts=c(1)
  for(j in 2:ncol(lirt)) lcuts=c(lcuts,max(lcuts)+ifelse(lirt[1,j]<lirt[2,j-1],0,1))
  ldups[[i]]=tapply(colnames(lirt),lcuts,list)
    }
    ldups=ldups2=unlist(ldups,rec=F)
  ndups=sapply(ldups,length)
    }

  if(any(ndups>1)){
    for(i in which(ndups>1)){
      tmp=dfeic[dfeic$eic%in%ldups[[i]],]
      tmp=tmp[tmp$eic%in%names(which(rowSums(table(tmp$eic,tmp$samp))>5)),]
      if(length(unique(tmp$eic))<2){ldups2[[i]]=as.list(ldups[[i]]);next}
      xnew=seq(floor(min(tmp[, whichrt] * 10))/10, ceiling(max(tmp[,whichrt] * 10))/10, eicParams$bw/2)
      amts=tapply(1:nrow(tmp),tmp$eic,function(x)
        .GRconvEIC(tmp[x,],whichrt=whichrt,bw=eicParams$nsmo1*eicParams$bw,delrt=eicParams$bw/2,xnew=xnew))
      apcs=sapply(amts,function(x) prcomp(x[[1]],center=F)$x[,1])
      lcl=cutree(hclust(as.dist(1-abs(cor(apcs))),method="average"),h=.1)
      ldups2[[i]]=tapply(names(lcl),lcl,c)
      # drt=do.call("rbind",tapply(tmp[,whichrt],tmp$eic,range))
      # ldups2[[i]]=.GRmergeInterval(drt)
    }
    ldups2=unlist(ldups2,rec=F)
  }
  ndups2=sapply(ldups2,length)
  names(ldups2)=sapply(ldups2,function(x) x[1])
  
  ###
  dfeic2=dfeic
#  print(ndups2)
  if(any(ndups2>1)) 
    for(ix in which(ndups2>1)){
    if(verbose) cat(" merging:",ldups2[[ix]],"\n",sep=" ")
       ## fix ineic
      l=which(dfeic2$eic%in%ldups2[[ix]])
      tmp=dfeic2[l,]
      neic=names(which.max(table(tmp$eic)))
      if("ineic"%in%names(tmp)){
        tmpsc=apply(as.matrix(tmp[,c("scan","samp","mz")]),1,paste,collapse=";;")
      tmp[,"ineic"]=tapply(tmp[,"ineic"],tmpsc,max)[tmpsc]
      neic=names(which.max(table(tmp$eic[tmp["ineic"]>0])))
      }
      tmp=tmp[which(!duplicated(tmp[,c("scan","samp","mz")])),]
      tmp$eic=neic
      dfeic2=dfeic2[-l,]
      dfeic2=rbind(dfeic2,tmp)
    }

  ######
  # remove duplicated and with na
  l=1:nrow(dfeic2)
  isdups=tapply(l,dfeic2$eic[l],function(x) sum(table(dfeic2$scan[x],dfeic2$samp[x])>1))
  for(ix in names(which(isdups>0))){
    lix=which(dfeic2$eic==ix)
    tmpx=dfeic2[lix,]
    scca=apply(tmpx[,c("scan","samp")],1,paste,collapse=";")
    idups=names(which(table(scca)>1))
    if(length(idups)>0){
      l=which(scca%in%idups)
      l=l[which(is.na(tmpx[l,"mz"]))]
      if(length(l)>0) dfeic2=dfeic2[-lix[l],]
    }
  }
  
  
  l=1:nrow(dfeic2)
  isdups=tapply(l,dfeic2$eic[l],function(x) sum(table(dfeic2$scan[x],dfeic2$samp[x])>1))
  
  ### if there is gap>eicParams$dPPM -> split!
  for(ix in names(which(isdups>0))){
    lix=which(dfeic2$eic==ix)
    tmpx=dfeic2[lix,]
    lspmz=.GRsplist(tmpx[,whichmz],d=eicParams$dPPM,ismass=T)
    if(length(lspmz)>1){
      spmz=sapply(lspmz,function(x) range(tmpx[x,whichmz],na.rm=T))
      spmz[,1]=(1-0.3*10^-6*eicParams$dPPM/spmz[,1])*spmz[,1]
      spmz[,2]=(1+0.3*10^-6*eicParams$dPPM/spmz[,2])*spmz[,2]
      for(ik in 1:nrow(spmz)){
        l=lix[which(dfeic2[lix,whichmz]>=spmz[ik,1] & dfeic2[lix,whichmz]<=spmz[ik,2])]
        dfeic2[l,]$eic=paste(ix,"-",ik,sep="")
      }
    }
  }
  dfeic2=dfeic2[order(dfeic2$eic,dfeic2$samp,dfeic2$scan,dfeic2[,whichmz]),]
  return(dfeic2)
}
