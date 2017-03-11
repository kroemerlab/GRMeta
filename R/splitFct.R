
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


################
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

################

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


################

.GRsplitBoth2<-function(imz,irt,dppm=10,dmz=NA,drt=1.1,typ="max"){
  l0=list(1:length(imz))
  doLoop=TRUE
  while(doLoop){
    l=unlist(lapply(l0,function(x) .GRsplistMZ(v = imz[x],iv = x,dppm = dppm,dmz = dmz,typ=typ)),recursive = F)
    l=unlist(lapply(l,function(x) .GRsplist(v = irt[x],iv = x,d=drt,ismass = F)),recursive = F)
    if(length(l)==length(l0)) doLoop=FALSE
    l0=l
  }
  # chk=t(sapply(l0,function(x) c(range(imz[x]),range(irt[x]))))
  # chk[order(chk[,1]),,drop=F]
  return(l0)
}

################
.GRsplistMZ<-function (v, iv = 1:length(v), dppm = NA,dmz=NA,typ="max"){
  lso = order(v)
  v2 = v[lso]
  lsup=NULL
  if(!is.na(dmz)){
    cdppm=10^6*dmz/v2[-length(v2)]
    if(!is.na(dppm))
      if(typ=="min") cdppm=cdppm %>% {.[.>dppm] = dppm; .} else cdppm=cdppm %>% {.[.<dppm] = dppm; .}
  } else  cdppm=rep(dppm,length(lso)-1)
  
  lsup=which(.GRcompdppm(v2,FALSE)>cdppm)
  
  lend = c(lsup, length(v))
  lst = c(1, lsup + 1)
  lapply(1:length(lst), function(x) iv[lso[lst[x]:lend[x]]])
}

