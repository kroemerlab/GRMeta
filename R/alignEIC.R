.GRcompAlignOneEIC<-function(ieic,lrefs,tabeic,eicfile=NULL,whichrt="rtcor",
                       eicParams){
  
  require(mgcv)

if(is.null(eicfile)) eicfile=paste(eicParams$dirEic,tabeic$GrpEic[which(tabeic$Id==ieic)],
                                   eicParams$addEic,".rda",sep="")

  load(eicfile)
  
  tmpeic=dfeic[dfeic$eic==ieic,]
  meic = .GRconvEIC(tmpeic, whichrt = whichrt, bw = eicParams$nsmo1 * eicParams$bw, delrt = eicParams$bw/2)
        
  m=meic[[1]]           
  xrt=meic[[2]]           
  
  rtref=codaref=rep(NA,length(lrefs));names(rtref)=names(codaref)=lrefs
  cmax=drt=matrix(NA,nrow=length(lrefs),ncol=ncol(m),dimnames=list(lrefs,colnames(m)))

  lenout=2^ceiling(log2(nrow(m)*2))
   for(iref in lrefs[lrefs%in%colnames(m)]){
    rec=apply(m,2,.GRcompVals,m[,iref],lenout,21)
  cmax[iref,]=apply(rec,2,max)
  drt[iref,]=(which.max(rec[,iref])-apply(rec,2,which.max))*median(diff(rt))
  rtref[iref]=weighted.median(xrt,m[,iref])
  codaref[iref]= .GRcodadw2(tmpeic$y[tmpeic$samp==iref])
}
  return(list(cmax=cmax,drt=drt,rtref=rtref,codaref=codaref))
  
}

################
.GRcompAlignGrpEIC<-function(grpeic,lrefs,tabeic=NULL,eicfile=NULL,whichrt="rtcor",eicParams){
  
  
  if(is.null(eicfile)) eicfile=paste(eicParams$dirEic,grpeic,
                                     eicParams$addEic,".rda",sep="")
  
  load(eicfile)
  if(is.null(tabeic)) leics=tabeic$Id[tabeic==grpeic] else leics=unique(dfeic$eic)
  leics=leics[leics%in%dfeic$eic]
  
  ares=list()
  for(ieic in leics){
  tmpeic=dfeic[dfeic$eic==ieic,]
  meic = .GRconvEIC(tmpeic, whichrt = whichrt, bw = eicParams$nsmo1 * eicParams$bw, delrt = eicParams$bw/2)
  
  m=meic[[1]]           
  xrt=meic[[2]]           
  
  rtref=codaref=rep(NA,length(lrefs));names(rtref)=names(codaref)=lrefs
  cmax=drt=matrix(NA,nrow=length(lrefs),ncol=ncol(m),dimnames=list(lrefs,colnames(m)))
  
  lenout=2^ceiling(log2(nrow(m)*2))
  for(iref in lrefs[lrefs%in%colnames(m)]){
    rec=apply(m,2,.GRcompVals,m[,iref],lenout,21)
    cmax[iref,]=apply(rec,2,max)
    drt[iref,]=(which.max(rec[,iref])-apply(rec,2,which.max))*median(diff(rt))
    rtref[iref]=weighted.median(xrt,m[,iref])
    codaref[iref]= .GRcodadw2(tmpeic$y[tmpeic$samp==iref])
  }
  ares[[ieic]]=list(cmax=cmax,drt=drt,rtref=rtref,codaref=codaref)
  }
  return(ares)
}


.GRcompAlignMGrpEIC<-function(tabeic,lrefs,whichrt="rtcor",eicParams,ncl=1){
                             
# 
#   if(is.null(names(lfiles))) names(lfiles)=1:length(lfiles)
#   lfiles=lfiles[file.exists(lfiles)]
#   if(length(lfiles)==0) return(NULL)
#   
  require(mgcv)
  
  leicgrp=unique(tabeic$GrpEic)
  leicgrp=leicgrp[file.exists(paste(eicParams$dirEic,leicgrp,eicParams$addEic,".rda",sep=""))]
  
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
  }
  d0=proc.time()[3]
  cat("Started at ",date(),sep="")
  
  acef=list()
  if(ncl==1){
    cat(" on 1 processor\n",sep="")
    acef=lapply(leicgrp,.GRcompAlignGrpEIC,
                          lrefs,tabeic,eicfile=NULL,whichrt=whichrt,eicParams)
  }
  if(ncl>1){
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfLibrary(GRMeta)
    sfLibrary(mgcv)
    #if(local) 
#    sfExport( ".GRcompAlignGrpEIC", local=TRUE )
    acef=sfClusterApplyLB(leicgrp,.GRcompAlignGrpEIC,
                          lrefs,tabeic,eicfile=NULL,whichrt=whichrt,eicParams)
    sfStop()
  }  
  acef=unlist(acef,recursive=FALSE)
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(leicgrp),1)," secs per file\n",sep="")
  return(invisible(acef))
  
  
}

# .GRcompAlignEICs

# 
# lsamps=unique(unlist(lapply(ares,function(x) colnames(x[[1]]))))
# lrefs=unique(unlist(lapply(ares,function(x) rownames(x[[1]]))))
# 
# adrt=lapply(ares,function(x){
#   x=x$drt[match(lrefs,rownames(x$drt)),match(lsamps,colnames(x$drt))]
#   dimnames(x)=list(lrefs,lsamps);x
#   })
# acmax=lapply(ares,function(x){
#   x=x$cmax[match(lrefs,rownames(x$cmax)),match(lsamps,colnames(x$cmax))]
#   dimnames(x)=list(lrefs,lsamps);x
# })
# 
# artref=sapply(ares,function(x) x$rtref[lrefs])
# acodaref=sapply(ares,function(x) x$codaref[lrefs])
# 
# newrts=(floor(min(artref,na.rm=T)/(2*eicParams$bw))-2):(ceiling(max(artref,na.rm=T)/(2*eicParams$bw))+2)
# newrts=newrts*eicParams$bw*2
# 
# ###################
# # for each sample
# allmpr=list()
# for(isamp in lsamps){
#   allmpr[[isamp]]=matrix(NA,nrow=length(newrts),ncol=length(lrefs),dimnames=list(newrts,lrefs))
# mrt=sapply(adrt,function(x) x[,isamp])
# cmax=sapply(acmax,function(x) x[,isamp])
# cmax[which(cmax<0.2)]=0
# 
# for(iref in lrefs){
# df=na.omit(data.frame(y=mrt[iref,],x=artref[iref,],w=cmax[iref,],cw=acodaref[iref,]))
# df$w[df$w<0.2]=0
# if(sum(df$w>0)>11){
#   m<-try(mgcv:::gam(y~s(x),data=df,weights=w),TRUE)
#   if("try-error"%in%class(m))  m<-try(mgcv:::gam(y~(x),data=df,weights=w),TRUE)
#   pr=mgcv:::predict.gam(m,data.frame(x=newrts))
#   #  pr[rts1<min(df$x,na.rm=T) | rts1>max(df$x,na.rm=T)]=NA
#   allmpr[[isamp]][,iref]=pr
# }
# }
# }
# 
# iref
# 
# 
# mpr=matrix(0,nrow=length(rts1),ncol=length(lsamp),dimnames=list(rts1,lsamp))
# for(isamp in lsamp[lsamp!=iref]){
#   df=na.omit(data.frame(y=msh[,isamp],x=refrt,w=mcor[,isamp]))
#   df$w[df$w<0.2]=0
#   if(sum(df$w>0)>11){
#     m<-try(mgcv:::gam(y~s(x),data=df,weights=w),TRUE)
#     if("try-error"%in%class(m))  m<-try(mgcv:::gam(y~(x),data=df,weights=w),TRUE)
#     pr=mgcv:::predict.gam(m,data.frame(x=rts1))*delrt
#     #  pr[rts1<min(df$x,na.rm=T) | rts1>max(df$x,na.rm=T)]=NA
#     mpr[,isamp]=pr
#   }
# }
# 
# 
# 
# 
#   leics=rownames(tousemat)
#   lsamp=colnames(tousemat)
#   mcor=msh=matrix(NA,ncol=length(lsamp),nrow=length(leics),dimnames=list(leics,lsamp))
#   refrt=rep(NA,length(leics));names(refrt)=leics
#   leics2=names(which(tousemat[,as.numeric(iref)]))
#   d0=proc.time()[3]
#   
#   if(serial) cat("Reference ",iref," (",length(leics2),")",sep="")  
#   for(i in 1:length(leicf)){
#     load(leicf[i])
#     if(serial) cat(".")
#     leic=as.character(unique(dfeic$Eic))
#     leic=leic[leic%in%leics2]
#     for(ix in leic){
#       tt=dfeic[dfeic$Eic==ix,]
#       vmint=filP$Mint[ix]
#       if(is.na(vmint)) vmint=median(filP$Mint)
#       re<-filtreScan(csc=as.vector(tt[,"y"]>vmint),idx=tt$Samp,alpha=filP$LSc,perc=filP$Perc)
#       tt$isin=re
#       
#       lrt=rts1[which(rts1>=(min(tt[,whichrt])-delrt) & rts1<=(max(tt[,whichrt])+delrt))]
#       if(max(lrt)<max(tt[,whichrt])) lrt=c(lrt,max(lrt)+(1:ceiling((max(tt[,whichrt])-max(lrt))/delrt))*delrt)
#       if(min(lrt)>min(tt[,whichrt])) lrt=c(min(lrt)-(ceiling((min(lrt)-min(tt[,whichrt]))/delrt):1)*delrt,lrt)
#       m=sapply(unique(tt$Samp), function(ik) approx(tt[,whichrt][tt$Samp==ik],tt$y[tt$Samp==ik],lrt)$y)
#       dimnames(m)=list(lrt,unique(tt$Samp))
#       m[is.na(m)]=0
#       
#       lenout=2^ceiling(log2(nrow(m)*2))
#       rec=apply(m,2,compVals,m[,iref],lenout,maxsh)
#       #  matplot(rec,typ="l",main=paste(iref,ix),col=mfiles$cols[as.numeric(colnames(rec))],lwd=2,lty=1)
#       mcor[ix,]=apply(rec,2,max)[as.character(lsamp)]
#       msh[ix,]=apply(rec,2,which.max)[as.character(lsamp)]-which.max(rec[,iref])
#       if(sum(tt$Samp==iref & tt$isin)==0) print(c(i,ix))
#       refrt[ix]=weighted.median(tt[,whichrt][tt$Samp==iref & tt$isin],tt$y[tt$Samp==iref & tt$isin])
#     }
#   }
#   d1=proc.time()[3]
#   if(serial) cat(":",round(d1-d0,1),"secs.\n Compute GAMs")
#   top=paste("Reference ",iref," (",length(leics2),"): ",round(d1-d0,1),"secs.",sep="")
#   #  save(file=paste("totot",iref,".rda",sep=""),refrt,msh,mcor,lsamp,mfiles,leicf,tousemat,maxsh,whichrt,filP,border,rts1,delrt)
#   
#   mpr=matrix(0,nrow=length(rts1),ncol=length(lsamp),dimnames=list(rts1,lsamp))
#   for(isamp in lsamp[lsamp!=iref]){
#     df=na.omit(data.frame(y=msh[,isamp],x=refrt,w=mcor[,isamp]))
#     df$w[df$w<0.2]=0
#     if(sum(df$w>0)>11){
#       m<-try(mgcv:::gam(y~s(x),data=df,weights=w),TRUE)
#       if("try-error"%in%class(m))  m<-try(mgcv:::gam(y~(x),data=df,weights=w),TRUE)
#       pr=mgcv:::predict.gam(m,data.frame(x=rts1))*delrt
#       #  pr[rts1<min(df$x,na.rm=T) | rts1>max(df$x,na.rm=T)]=NA
#       mpr[,isamp]=pr
#     }
#   }
#   ### smooth border
#   #  save(file=paste("totot",iref,".rda",sep=""),refrt,mpr)
#   
#   lb1=which(rts1>min(refrt,na.rm=T) & rts1<=(min(refrt,na.rm=T)+border+2*delrt))
#   l2rep=1:lb1[(floor(length(lb1)/2))]
#   mpr[l2rep,]=matrix(rep(apply(mpr[lb1,],2,median),each=length(l2rep)),ncol=length(lsamp))
#   lb1=which(rts1<max(refrt,na.rm=T) & rts1>=(max(refrt,na.rm=T)-border-2*delrt))
#   l2rep=nrow(mpr):rev(lb1)[(floor(length(lb1)/2))]
#   mpr[l2rep,]=matrix(rep(apply(mpr[lb1,],2,median),each=length(l2rep)),ncol=length(lsamp))
#   if(serial) cat("\n")
#   
#   return(list(sht=msh,pr=mpr,cor=mcor,rts=rts1,ref=iref,top=top))
# }
# 
# 
# #######
# ## comp alignment based on EIC
# compAlignMEIC<-function(tabeic,lrefs,whichrt="rtcor",maxsh=50,params,border=0.2,ncl=1){
#   
#   require(mgcv)
# 
#   d0=proc.time()[3]
#   
#   if(missing(lrefs)) stop('Sid of the ref samples')
#   
#   
#   if(ncl>1){
#     if(ncl>parallel:::detectCores()) ncl=min(ncl,parallel:::detectCores()-1,length(lrefs))
#     if(sfIsRunning()) sfStop()
#     sfInit(parallel=TRUE, cpus=ncl, type="SOCK")
#     sfExport("filtreScan","compVals", "compAlignEIC","compAlignEICCl", local=TRUE )
#     allr=sfClusterApplyLB(1:length(lrefs),compAlignEICCl,lrefs,mfiles,leicf,tousemat,maxsh,whichrt,filP,border)
#     names(allr)=lrefs
#     sfStop()
#     allr=allr[sapply(allr,length)>0]
#     cat(sapply(allr,function(x) x$top),sep="\n")
#   }
#   
#   if(ncl==1){
#     allr=list()
#     for(iref in lrefs){
#       resp<-try(compAlignEIC(iref=iref,mfiles,leicf,tousemat,maxsh=maxsh,whichrt=whichrt,filP=filP,border=border,serial=TRUE),TRUE)
#       if(class(resp)!="try-error") allr[[iref]]<-resp
#     }
#   }
#   d1=proc.time()[3]
#   
#   cat("\n Took",round(d1-d0,1),"secs. to process ",length(allr)," refs\n")
#   
#   return(list(allres=allr,whichrt=whichrt,border=border))
# }