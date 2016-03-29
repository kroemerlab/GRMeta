.inGRconvMS<-function(lfiles,lsids,
                      ldatas=c(MZ="MZap",Area="Area.scan",AreaSm="Area.scansm",
                               Height="HE.ap",RT="RT.ap",RTsm="RT.apsm",Coda="Coda",NCons="NCons")){
  alleics=alldata=list()
  for(ifi in 1:length(lfiles)){
    system.time(load(lfiles[ifi]))
    rm("dfeic")
    alldat=lapply(alldat,function(x) x[match(lsids,rownames(x)),,drop=F])
    alldat=alldat[names(alldat)%in%ldatas]
    for(i in names(alldat)) if(is.null(alldata[[i]])) alldata[[i]]=alldat[[i]] else alldata[[i]]=cbind(alldata[[i]],alldat[[i]])
    alleics[[ifi]]=eicstats
  }
  alldata=alldata[!sapply(alldata,is.null)]
  eicsdat=do.call("rbind",alleics)
  rm(list=c("alleics",'alldat','eicstats'))
  
  leics=eicsdat$PkId
  alldata=lapply(alldata,function(x) x[,match(leics,colnames(x)),drop=F])
  names(alldata)=names(ldatas[match(names(alldata),ldatas)])
  for(k in names(alldata)) dimnames(alldata[[k]])=list(lsids,leics)
  return(list(dat=alldata,eic=eicsdat))
}

