.GRcompAlignOneEIC<-function(ieic,lrefs,tabeic,eicfile=NULL,whichrt="rtcor",
                       eicParams){
  
  require(mgcv)

if(is.null(eicfile)) eicfile=paste(eicParams$dirEic,tabeic$GrpEic[which(tabeic$Id==ieic)],
                                   eicParams$addEic,".rda",sep="")

  load(eicfile)
  tmpeic=dfeic[dfeic$eic==ieic,]
  
  rtref=codaref=rep(NA,length(lrefs));names(rtref)=names(codaref)=lrefs
  cmax=drt=matrix(NA,nrow=length(lrefs),ncol=ncol(m),dimnames=list(lrefs,colnames(m)))
  if(any(!lrefs%in%dfeic$samp)) return(list(cmax=cmax,drt=drt,rtref=rtref,codaref=codaref))
  
  meic = .GRconvEIC(tmpeic, whichrt = whichrt, bw = eicParams$nsmo1 * eicParams$bw, delrt = eicParams$bw/2)
        
  m=meic[[1]]           
  xrt=meic[[2]]           
  
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
.GRcompAlignGrpEIC<-function(grpeic,lrefs,drtmax=0.3,tabeic=NULL,eicfile=NULL,whichrt="rtcor",eicParams,verbose=FALSE){
  
  require(limma)
  if(is.null(eicfile)) eicfile=paste(eicParams$dirEic,grpeic,
                                     eicParams$addEic,".rda",sep="")
  if(verbose) cat("* ",eicfile,": ")
  load(eicfile)
  if(!is.null(tabeic)) leics=tabeic$Id[tabeic==grpeic] else leics=unique(dfeic$eic)
  leics=leics[leics%in%dfeic$eic]
  
  ares=list()
  for(ieic in leics){
    if(verbose) cat(ieic," ",sep=" ")
    tmpeic=dfeic[dfeic$eic==ieic,]
  meic = .GRconvEIC(tmpeic, whichrt = whichrt, bw = eicParams$nsmo1 * eicParams$bw, delrt = eicParams$bw/2)
  
  m=meic[[1]]           
  xrt=meic[[2]]           
  
  rtref=codaref=rep(NA,length(lrefs));names(rtref)=names(codaref)=lrefs
  cmax=drt=matrix(NA,nrow=length(lrefs),ncol=ncol(m),dimnames=list(lrefs,colnames(m)))
  
  nmax=ceiling(drtmax/median(diff(xrt))+3)
  if(nmax%%2==0) nmax=nmax+1
  lenout=2^ceiling(log2(nrow(m)*2))
  for(iref in lrefs[lrefs%in%colnames(m)]){
    rec=apply(m,2,.GRcompVals,m[,iref],lenout,nmax)
    cmax[iref,]=apply(rec,2,max)
    drt[iref,]= (apply(rec,2,which.max)-which.max(rec[,iref]))*median(diff(xrt))
    rtref[iref]=weighted.median(xrt,m[,iref])
    codaref[iref]= .GRcodadw2(tmpeic$y[tmpeic$samp==iref])
  }
  ares[[ieic]]=list(cmax=cmax,drt=drt,rtref=rtref,codaref=codaref)
  }
  if(verbose) cat("\n")
  
  return(ares)
}

################
.GRcorGrpEIC<-function(lgrpeic,shmat,tabeic=NULL,byGrp=FALSE,dosave=FALSE,
                       whichrt="rtcor",newrt="rtcor2",eicParams,verbose=FALSE,retres=TRUE,ncl=1){
  
  .inGRcorGrpEIC<-function(grpeic,shmat,tabeic=NULL,byGrp=FALSE,dosave=FALSE,
                           whichrt="rtcor",newrt="rtcor2",eicParams,verbose=TRUE,retres=TRUE){
    
  eicfile=paste(eicParams$dirEic,grpeic,eicParams$addEic,".rda",sep="")
  if(verbose) cat("* ",eicfile)
  load(eicfile)
  if(!is.null(tabeic)) leics=tabeic$Id[tabeic$GrpEic==grpeic] else leics=unique(dfeic$eic)
  leics=leics[leics%in%dfeic$eic]
  
  eicinf=attr(dfeic,"eic")
  eicinf=eicinf[eicinf$Id%in%leics,]
  rts=as.numeric(rownames(shmat))
  addrt=apply(shmat,2,function(x) approx(rts,x,eicinf$rtap)$y)
  rownames(addrt)=eicinf$Id
  addrtv=as.vector(addrt)
  names(addrtv)=paste(rep(rownames(addrt),ncol(addrt)),rep(colnames(addrt),each=nrow(addrt)),sep=";;;")
  
  addnewrt=addrtv[paste(dfeic$eic,dfeic$samp,sep=";;;")]
  if(any(is.na(addnewrt))){
    if(verbose) cat(" some missing!")
    addnewrt[is.na(addnewrt)]=0
  }
  addnewrt=dfeic[,whichrt]+addnewrt
  dfeic[,newrt]=addnewrt
  
  if(dosave) save(file=eicfile,dfeic,eicinfos,sampinfos)
  
  if(verbose) cat("\n")
  if(retres) invisible(dfeic) else return(nrow(eicinf))
  }
  
  ###################
  if(ncl!=1){
    require("snowfall")
    ncl=max(1,min(ncl,parallel:::detectCores()))
  }
  
  if(!is.null(tabeic) & is.null(lgrpeic))  lgrpeic=unique(tabeic$GrpEic)
  lfiles=paste(ifelse(is.null(eicParams$dirEic),"./",eicParams$dirEic),
               lgrpeic,
               ifelse(is.null(eicParams$addEic),"./",eicParams$addEic),".rda",sep="")
  lgrpeic=unique(lgrpeic[file.exists(lfiles)])
  if(!is.null(tabeic)) lgrpeic=unique(tabeic$GrpEic[tabeic$GrpEic%in%lgrpeic])
  
  d0=proc.time()[3]
  cat("Started at ",date()," on ",ncl," processors\n",sep="")
  if(ncl==1){
    allr=lapply(lgrpeic,.inGRcorGrpEIC,shmat=shmat,tabeic=tabeic,byGrp=byGrp,dosave=dosave,whichrt=whichrt,newrt=newrt,eicParams,
                verbose=verbose,retres=retres)
   if(retres) totop=sum(sapply(allr,function(x) length(unique(dfeic$eic)))) else totop=sum(unlist(allr))
  }
  if(ncl>1){
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfExport( ".inGRcorGrpEIC", local=TRUE )
    sfLibrary(GRMeta)
    sfLibrary(limma)
    allr=sfClusterApplyLB(lgrpeic,.inGRcorGrpEIC,shmat=shmat,tabeic=tabeic,byGrp=byGrp,dosave=TRUE,
                          whichrt=whichrt,newrt=newrt,eicParams,verbose=FALSE,retres=FALSE)
    totop=sum(unlist(allr))
    sfStop()
  }
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs to process ",length(lgrpeic)," EIC groups\n",sep="")
  cat(" * ",totop," EICs in all, ",round(totop/length(lgrpeic),3)," on av.\n")
  return(invisible(allr))
}

################
.GRcorGrpEICo<-function(grpeic,shmat,tabeic=NULL,eicfile=NULL,byGrp=FALSE,dosave=FALSE,
                       whichrt="rtcor",newrt="rtcor2",eicParams,verbose=FALSE,ncl=1){
  
  if(is.null(eicfile)) eicfile=paste(eicParams$dirEic,grpeic,
                                     eicParams$addEic,".rda",sep="")
  if(verbose) cat("* ",eicfile)
  load(eicfile)
  if(!is.null(tabeic)) leics=tabeic$Id[tabeic==grpeic] else leics=unique(dfeic$eic)
  leics=leics[leics%in%dfeic$eic]
  
  eicinf=attr(dfeic,"eic")
  eicinf=eicinf[eicinf$Id%in%leics,]
  rts=as.numeric(rownames(shmat))
  addrt=apply(shmat,2,function(x) approx(rts,x,eicinf$rtap)$y)
  rownames(addrt)=eicinf$Id
  addrtv=as.vector(addrt)
  names(addrtv)=paste(rep(rownames(addrt),ncol(addrt)),rep(colnames(addrt),each=nrow(addrt)),sep=";;;")
  
  addnewrt=addrtv[paste(dfeic$eic,dfeic$samp,sep=";;;")]
  if(any(is.na(addnewrt))){
    if(verbose) cat(" some missing!")
  addnewrt[is.na(addnewrt)]=0
  }
  addnewrt=dfeic[,whichrt]+addnewrt
  dfeic[,newrt]=addnewrt
  
  if(dosave) save(file=eicfile,dfeic,eicinfos,sampinfos)

  if(verbose) cat("\n")
  

  return(invisible(dfeic))
}

.GRcompAlignMGrpEIC<-function(tabeic,lrefs,drtmax=0.3,whichrt="rtcor",eicParams,ncl=1){
                             
# 
#   if(is.null(names(lfiles))) names(lfiles)=1:length(lfiles)
#   lfiles=lfiles[file.exists(lfiles)]
#   if(length(lfiles)==0) return(NULL)
#   
  require(limma)
  
  leicgrp=unique(tabeic$GrpEic)
  leicgrp=leicgrp[file.exists(paste(eicParams$dirEic,leicgrp,eicParams$addEic,".rda",sep=""))]
  print(leicgrp)
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
                          lrefs,tabeic,drtmax=drtmax,eicfile=NULL,whichrt=whichrt,eicParams,verbose=TRUE)
  }
  if(ncl>1){
    cat(" on ",ncl," processors\n",sep="")
    sfInit(parallel=TRUE, cpus=ncl, type='SOCK')
    sfLibrary(GRMeta)
    sfLibrary(limma)
    #if(local) 
#    sfExport( ".GRcompAlignGrpEIC", local=TRUE )
    acef=sfClusterApplyLB(leicgrp,.GRcompAlignGrpEIC,
                          lrefs,tabeic,drtmax=drtmax,eicfile=NULL,whichrt=whichrt,eicParams)
    sfStop()
  }  
  acef=unlist(acef,recursive=FALSE)
  d1=proc.time()[3]
  cat("\nCompleted at ",date()," -> took ",round(d1-d0,1)," secs - ",round((d1-d0)/length(leicgrp),1)," secs per file\n",sep="")
  return(invisible(acef))
  
  
}

.GRcompAlignEICs<-function(ares,newrts=NULL,lsamps=NULL,lrefs=NULL,rtlim=NULL,drtmax=0.5){


if(is.null(lsamps)) lsamps=unique(unlist(lapply(ares,function(x) colnames(x[[1]]))))
if(is.null(lrefs)) lrefs=unique(unlist(lapply(ares,function(x) rownames(x[[1]]))))

adrt=lapply(ares,function(x){
  x=x$drt[match(lrefs,rownames(x$drt)),match(lsamps,colnames(x$drt))]
  dimnames(x)=list(lrefs,lsamps);x
  })
acmax=lapply(ares,function(x){
  x=x$cmax[match(lrefs,rownames(x$cmax)),match(lsamps,colnames(x$cmax))]
  dimnames(x)=list(lrefs,lsamps);x
})

artref=sapply(ares,function(x) x$rtref[lrefs])
acodaref=sapply(ares,function(x) x$codaref[lrefs])

if(is.null(newrts)){
  toadd=max(2,ceiling(drtmax/eicParams$bw))
  newrts=(floor(min(artref,na.rm=T)/(2*eicParams$bw))-toadd):(ceiling(max(artref,na.rm=T)/(2*eicParams$bw))+toadd)
newrts=newrts*eicParams$bw*2
}

doborderu=doborderb=FALSE
if(!is.null(rtlim) & any(newrts<rtlim[1])){
  lbe=which(newrts<rtlim[1])
  lbe2=which(newrts>=rtlim[1])[1:3]
  doborderb=TRUE
}
if(!is.null(rtlim) & any(newrts>rtlim[2])){
  lup=which(newrts>rtlim[2])
  lup2=rev(which(newrts<=rtlim[2]))[1:3]
  doborderu=TRUE
}


###################
# for each sample
allmpr=list()
for(isamp in lsamps){
  cat(".")
  allmpr[[isamp]]=matrix(NA,nrow=length(newrts),ncol=length(lrefs),dimnames=list(newrts,lrefs))
mrt=sapply(adrt,function(x) x[,isamp])
cmax=sapply(acmax,function(x) x[,isamp])
cmax[which(cmax<0.2)]=0

for(iref in lrefs){
df=na.omit(data.frame(y=mrt[iref,],x=artref[iref,],w=cmax[iref,],cw=acodaref[iref,]))
df$w[df$w<0.2]=0
if(sum(df$w>0)>11){
  m<-try(mgcv:::gam(y~s(x),data=df,weights=w),TRUE)
  if("try-error"%in%class(m))  m<-try(mgcv:::gam(y~(x),data=df,weights=w),TRUE)
  pr=mgcv:::predict.gam(m,data.frame(x=newrts))
  if(any(abs(pr)>drtmax)) pr[which(abs(pr)>drtmax)]=drtmax*sign(pr[which(abs(pr)>drtmax)])
  if(doborderb)  pr[lbe]=median(pr[lbe2])
  if(doborderu)  pr[lup]=median(pr[lup2])
  allmpr[[isamp]][,iref]=pr
}
}
}
cat("\n")

shmat=sapply(allmpr,function(x){
  pr=apply(x,1,median)
  if(any(abs(pr)>drtmax)) pr[which(abs(pr)>drtmax)]=drtmax*sign(pr[which(abs(pr)>drtmax)])
  if(doborderb)  pr[lbe]=median(pr[lbe2])
  if(doborderu)  pr[lup]=median(pr[lup2])
  pr=mgcv:::gam(pr~s(newrts))$fit
  if(doborderb)  pr[lbe]=pr[lbe2[which.min(lbe2)]]
  if(doborderu)  pr[lup]=pr[lup2[which.max(lup2)]]
  pr
})


rownames(shmat)=newrts

return(list(ShMat=shmat,RT=newrts,AllPr=allmpr,Ref=lrefs,RTref=artref))

}
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
