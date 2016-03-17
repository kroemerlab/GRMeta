partitionPeaks<-function(cdat,err.ppm=10,err.rt=0.1,metaboSet=list(),fac=4,whichrt="rt",whichmz="mz",whichsamp="samp"){
  
  #### Rough splitting unlikely neighbours: set to err.ppm*fac and err.rt*fac*2
  system.time(lsp<-.GRgenpkl(cdat[,whichmz],cdat[,whichrt],err.ppm=err.ppm*fac,err.rt=err.rt*fac*2))
  cat("Loose split: num. peaks:",nrow(cdat)," -> num. clusters: ",length(lsp),"\n",sep="")

  ##### Density based split on both rt and ppm
  system.time(lsp1<-lapply(lsp,function(x) 
    .GRsplitBoth(cdat[x,whichmz],cdat[x,whichrt],x,bwrt=err.rt,bwppm=err.ppm)))
  lsp1=unlist(lsp1,FALSE)
  lsp1=lapply(lsp1,function(x) x[order(cdat[x,whichsamp],cdat[x,whichsamp])])
  
  mrt=round(sapply(lsp1,function(x) median(cdat[x,whichrt])),1)
  mrz=sapply(lsp1,function(x) median(cdat[x,whichmz]))
  ana=sprintf("%.4f@%.2f",mrz,mrt)
  if(any(table(ana)>1)){
    l=which(ana%in%names(which(table(ana)>1)))
    ana[l]=sprintf("%.5f@%.3f",mrz,mrt)[l]
  }
  names(lsp1)=ana

  lso=order(mrt,mrz)
  lsp1=lsp1[lso]
  cat("Density based split: num. peaks=",nrow(cdat)," -> num. clusters:",length(lsp1),"\n",sep="")
  if(is.null(metaboSet))  return(lsp1)
  
  ##################
  ## conv2 metaboset
  # metaboSet=list(assay="cef",Meta=NULL,File=NULL)
  ## do annot
  mrt=sapply(lsp1,function(x) median(cdat[x,whichrt]))
  mrz=sapply(lsp1,function(x) median(cdat[x,whichmz]))
  ana=sprintf("%.4f@%.2f-%s",mrz,mrt,metaboSet$Method)
  if(any(table(ana)>1)){
    l=which(ana%in%names(which(table(ana)>1)))
    ana[l]=sprintf("%.5f@%.3f-%s",mrz,mrt,metaboSet$Method)[l]
  }
  names(lsp1)=ana
 
   Annot=data.frame(Analyte=ana,MetName=NA,IsSTD=FALSE,LevelAnnot=4,Method=metaboSet$Method,
                    MZ=mrz,RT=mrt,stringsAsFactors=FALSE)
  rownames(Annot)=ana
  lsids=unique(cdat[,whichsamp])
  ## do meta/files
  Meta=File=data.frame(Sid=lsids,sType="Sa",stringsAsFactors=F)
  Meta$InjOrder=1:nrow(Meta)
  File$Date=chron(1:nrow(Meta))
  rownames(Meta)=rownames(File)=lsids
  if(!is.null(metaboSet$Meta))
    if(all(lsids%in%metaboSet$Meta$Sid)) Meta=metaboSet$Meta[lsids,]
  if(!is.null(metaboSet$File))
    if(all(lsids%in%metaboSet$File$Sid)) File=metaboSet$File[lsids,]
  ## 
  mrt=sapply(lsp1,function(x) tapply(cdat[x,whichrt],cdat[x,whichsamp],median,na.rm=T)[lsids])
  mmz=sapply(lsp1,function(x) tapply(cdat[x,whichmz],cdat[x,whichsamp],median,na.rm=T)[lsids])
  dimnames(mrt)=dimnames(mmz)=list(lsids,names(lsp1))
  Data=list(RT=mrt,MZ=mmz)

  lvars=names(cdat)[sapply(cdat, class)%in%c("numeric","integer")]
  lvars=lvars[!lvars%in%c(whichrt,whichmz)]
  for(i in lvars){
    m=sapply(lsp1,function(x) tapply(cdat[x,i],cdat[x,whichsamp],sum,na.rm=T)[lsids])
    dimnames(m)=list(lsids,names(lsp1))
    Data[[i]]=m
  }

  allmat=list(Method=metaboSet$Method,Sid=lsids,Analyte=Annot$Analyte,Annot=Annot,Meta=Meta,File=File,Data=Data)
  class(allmat)=append(class(allmat),"metaboSet")
  invisible(allmat)
  
}


.GRgenpkl<-function(mass,rtime,err.ppm=10,err.rt=.1,err.mz=NULL){
  
  
  doone<-function(x,lsp,err.rt,err.ppm,ismass){
    np=sapply(lsp,length)
    if(ismass) dm=sapply(lsp,function(i) max(.GRcompdppm(x[i])))>err.ppm
    if(!ismass) dm=sapply(lsp,function(i) max(x[i])-min(x[i]))>err.rt
    l=which(dm & np>1)
    if(length(l)>0){
      d=ifelse(ismass,err.ppm,err.rt)
      lsp1=lapply(lsp[l],function(y) .GRsplist(x[y],y,ismass=ismass,d=d))
      if(any(sapply(lsp1,is.list))) lsp1=unlist(lsp1,FALSE)
      lsp=c(lsp[-l],lsp1)
    }
    lsp
  }
  
  clsp=list()
  if(is.null(err.mz)) lsp=.GRsplist(mass,ismass=TRUE,d=err.ppm)
  if(!is.null(err.mz)) lsp=.GRsplist(mass,ismass=FALSE,d=err.mz)
  i=0
  while(length(lsp)!=length(clsp)){
    clsp=lsp
    if((i%%2==1)){
      if(is.null(err.mz)) lsp=doone(mass,lsp,err.rt,err.ppm,ismass=TRUE)
      if(!is.null(err.mz)) lsp=doone(mass,lsp,err.mz,err.ppm,ismass=FALSE)
    }
    if((i%%2!=1)) lsp=doone(rtime,lsp,err.rt,err.ppm,ismass=FALSE)
    i=i+1
  }
  
  art=round(sapply(lsp,function(i) median(rtime[i])),2)
  amz=sapply(lsp,function(i) median(mass[i]))
  lsp=lsp[order(art,amz)]
  return(lsp)
  
}

.GRcompdppm<-function(v,sort=T){
  if(length(v)==1) return(0)
  if(sort) v=sort(v)
  drt=diff(v)/v[-length(v)]
  10^6*drt
}

.GRsplist<-function(v,iv=1:length(v),d=0.2,ismass=FALSE){
  lso=order(v)
  v2=v[lso]  
  if(ismass) drt=.GRcompdppm(v2) else drt=diff(v2)
  lend=c(which(drt>d),length(v))
  lst=c(1,which(drt>d)+1)
  lapply(1:length(lst),function(x) iv[lso[lst[x]:lend[x]]])
}


.GRsplitBoth<-function(mass,rtime,indices=1:length(mass),bwrt=0.1,bwppm=10){
  
  if(length(mass)==1) return(list(indices))
  
  ires=list(1:length(mass))
  cres=list()
  while(length(ires)>0){
    l=ires[[1]]
    resppm=.GRsplitOne(mass[l],bwppm,TRUE)
    resrt=.GRsplitOne(rtime[l],bwrt,FALSE)
    newcl=as.numeric(factor(max(resppm$cl)*resppm$cl+resrt$cl))
    res=tapply(l,newcl,c)
    if(length(res)==1) cres=c(cres,lapply(res,function(ix) indices[ix]))
    if(length(res)>1) ires=c(ires,res)
    ires=ires[-1]
  }
  return(cres)
}

.GRsplitOne<-function(v,bw,ismass=TRUE){
  
  medv=median(v)
  if(ismass) v=(v-medv)*10^6/medv
  
  nbins=2^ceiling(log2(round(diff(range(v))*10/bw)))
  nbins=max(512,nbins)
  dens=density(v,bw=bw,n=nbins)
  dens$nbins=nbins
  
  lpks=.GRfindturnpoint(dens$y)
  lpks$pks=lpks$pks[(dens$y[lpks$pks]/max(dens$y[lpks$pks]))>10^-6]
  irts=dens$x[lpks$pks]
  dirt=matrix(sapply(irts,function(irt) v-irt ),ncol=length(irts))
  newcl=apply(dirt^2,1,which.min)
  res=tapply(1:length(v),newcl,c)
  if(!is.list(res)) res=as.list(res)
  ### renorm
  if(ismass){
    irts=medv+irts*medv/10^6
    dens$x=medv+dens$x*medv/10^6
  }
  
  return(list(sp=res,cl=newcl,pk=irts,dens=dens))
}




.GRfindturnpoint<-function(y){
  peaks2 <- function(x, ties.method) {
    z <- embed(rev(as.vector(c(-Inf, x, -Inf))), dim = 3)
    z <- z[rev(seq(nrow(z))), ]
    v <- max.col(z, ties.method = ties.method) == 2
    v
  }
  msExtrema <- function(x) {
    l <- length(x)
    index1 <- peaks2(x, ties.method = "first")
    index2 <- peaks2(-x, ties.method = "last")
    index.max <- index1 & !index2
    index.min <- index2 & !index1
    list(index.max = index.max, index.min = index.min)
  }
  y <- y[!is.na(y)]
  if (length(unique(y)) == 1) {
    pks <- round(length(y)/2)
    vlys <- c(1, length(y))
    x <- new("list")
    x$pks <- pks
    x$vlys <- vlys
    return(x)
  }
  b <- msExtrema(y)
  pks <- which(b$index.max)
  vlys <- which(b$index.min)
  if (pks[1] != 1) 
    vlys <- c(1, vlys)
  if (pks[length(pks)] != length(y)) 
    vlys <- c(vlys, length(y))
  if (length(pks) == 1) 
    vlys <- c(1, length(y))
  x <- new("list")
  x$pks <- pks
  x$vlys <- vlys
  return(x)
}
