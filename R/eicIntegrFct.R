

######################################################################

.GRcompSeg<-function(cm,cmz,params,ivMint){
  
  pc=prcomp(cm,center=FALSE)
  ncps=max(min(ncol(pc$x),4),2)
  mm=pc$x[,1:ncps,drop=F]
  mmr=pc$r[,1:ncps,drop=F]
  mmz=(cmz%*%(mmr))
  
  
  #########
  
  v=1:nrow(mmz)
  lresid=apply(mmz,2,function(x) sum(residuals(lm(x~v-1))))
  lsi=sign(lresid)
  lsi[1]=sign(mm[which.max(abs(mm[,1])),1])
  lsi[lsi==0]=1
  mmr=sweep(mmr,2,lsi,"*")
  mm=sweep(mm,2,lsi,"*")
  mmz=sweep(mmz,2,lsi,"*")
  mnoise=colSums(mmr*ivMint)
  
  #####
  
  ipks=icmp=rep(0,nrow(mm))
  for(ipc in 1:2){
    v=(mm[,ipc])>abs(mnoise[ipc])
    if(any(v)){
      lst=which(diff(v)==1)-1
      if(length(lst)==0) lst=1
      len=which(diff(v)== -1)+1
      if(v[length(v)]) len=c(len,length(v))
      if(len[1]<lst[1]) lst=c(1,lst)
      for(i in 1:length(lst)){
        il=lst[i]:len[i]
        cpk=ipks[il]
        if(all(cpk==0)){
          ipks[il]=max(ipks)+1
          icmp[il]=ipc
        }
        if(length(unique(cpk[cpk!=0]))==1){
          ipks[il]=max(cpk)
          icmp[il]=max(icmp[il])
        }
        #if(!all(cpk==0) & !(length(unique(cpk[cpk!=0]))==1)) print("Error")
        
      }
    }
  }
  ### get apex
  ispk=rep(FALSE,length(ipks))
  lpks=names(which(table(ipks[ipks!=0])>(params$nspan)))
  if(length(lpks)>0){
    for(i in as.numeric(lpks)){
      lipk=range(which(ipks==i))
      iicmp=icmp[lipk[1]+1]
      lipk[1]=max(1,lipk[1]-2)
      lipk[2]=min(length(ispk),lipk[2]+2)
      lipk=lipk[1]:lipk[2]
      vm=mm[lipk,iicmp]
      vmz=mmz[lipk,iicmp]
      respk=.GRmsPeakSimple3(vm,noise.local=vmz,span=params$nspan,snr.thresh=1,compCl=F)
      if(nrow(respk)>0 & any(respk[,"aboveSN"]==1)){
        respk=respk[respk[,"aboveSN"]==1,c(1:4),drop=F]
        
        #### Group less< 12
        while(any(diff(respk[,1])<=(params$nspan*2+1))){
          lcl=cutree(hclust(dist(respk[,1,drop=F])),h=params$nspan*2+1)
          respk=respk[tapply(1:nrow(respk),lcl,function(x) x[which.max(respk[x,2])]),,drop=F]
        }
        
        #### Group less>=12 and 60 scans if 2 points and mininum<0.8 than the highest
        cond=((diff(respk[,1])>(params$nspan*2+1) & diff(respk[,1])<=(params$nspan*6+1)))
        chgt=T
        if(any(cond) & chgt){
          ch=FALSE
          dval=alp=cbind(respk[which(cond),1],respk[which(cond)+1,1])
          for(j in 1:nrow(dval)) dval[j,]=range(vm[alp[j,1]:alp[j,2]])
          drat=(dval[,2]-dval[,1])/dval[,2]
          #         print(drat)
          if(any(drat<0.2)){
            chgt=TRUE
            alp=alp[which.min(drat),]
            alp=alp[which.min(vm[alp])]
            respk=respk[respk[,1]!=alp,,drop=F]
            cond=((diff(respk[,1])>(params$nspan*2+1) & diff(respk[,1])<=(params$nspan*6+1)))
          }
        }
        
        ####
        ispk[lipk[respk[,1]]]=TRUE
      } ## if 
      
      
    } ## for i 
  } ## length lpks>0
  
  
  res=cbind(Segment=ipks,PCnum=icmp,Peak=ispk)
  return(cbind(res,mm[,1:ncps,drop=F]))
}

##########################################################

.GRcompSegOne<-function(cm,cmz,params,ivMint){
  
  ipks=icmp=rep(0,nrow(cm))
  v=(cm[,1]>cmz[,1])
  if(any(v)){
    lst=which(diff(v)==1)-1
    if(length(lst)==0) lst=1
    len=which(diff(v)== -1)+1
    if(v[length(v)]) len=c(len,length(v))
    if(len[1]<lst[1]) lst=c(1,lst)
    for(i in 1:length(lst)){
      il=lst[i]:len[i]
      cpk=ipks[il]
      if(all(cpk==0)){
        ipks[il]=max(ipks)+1
        icmp[il]=1
      }
      if(length(unique(cpk[cpk!=0]))==1){
        ipks[il]=max(cpk)
        icmp[il]=max(icmp[il])
      }
      #if(!all(cpk==0) & !(length(unique(cpk[cpk!=0]))==1)) print("Error")
      
    }
  }
  
  
  ### get apex
  ispk=rep(FALSE,length(ipks))
  lpks=names(which(table(ipks[ipks!=0])>(params$nspan)))
  if(length(lpks)>0){
    for(i in as.numeric(lpks)){
      lipk=range(which(ipks==i))
      iicmp=icmp[lipk[1]+1]
      lipk[1]=max(1,lipk[1]-2)
      lipk[2]=min(length(ispk),lipk[2]+2)
      lipk=lipk[1]:lipk[2]
      vm=cm[lipk,iicmp]
      vmz=cmz[lipk,iicmp]
      respk=.GRmsPeakSimple3(vm,noise.local=vmz,span=params$nspan,snr.thresh=1,compCl=F)
      if(nrow(respk)>0 & any(respk[,"aboveSN"]==1)){
        respk=respk[respk[,"aboveSN"]==1,c(1:4),drop=F]
        
        #### Group less< 12
        while(any(diff(respk[,1])<=(params$nspan*2+1))){
          lcl=cutree(hclust(dist(respk[,1,drop=F])),h=params$nspan*2+1)
          respk=respk[tapply(1:nrow(respk),lcl,function(x) x[which.max(respk[x,2])]),,drop=F]
        }
        
        #### Group less>=12 and 60 scans if 2 points and mininum<0.8 than the highest
        cond=((diff(respk[,1])>(params$nspan*2+1) & diff(respk[,1])<=(params$nspan*6+1)))
        chgt=T
        if(any(cond) & chgt){
          ch=FALSE
          dval=alp=cbind(respk[which(cond),1],respk[which(cond)+1,1])
          for(j in 1:nrow(dval)) dval[j,]=range(vm[alp[j,1]:alp[j,2]])
          drat=(dval[,2]-dval[,1])/dval[,2]
          #         print(drat)
          if(any(drat<0.2)){
            chgt=TRUE
            alp=alp[which.min(drat),]
            alp=alp[which.min(vm[alp])]
            respk=respk[respk[,1]!=alp,,drop=F]
            cond=((diff(respk[,1])>(params$nspan*2+1) & diff(respk[,1])<=(params$nspan*6+1)))
          }
        }
        
        ####
        ispk[lipk[respk[,1]]]=TRUE
      } ## if 
      
      
    } ## for i 
  } ## length lpks>0
  
  
  res=cbind(Segment=ipks,PCnum=icmp,Peak=ispk)
  v=cm[,1,drop=F]-cmz[,1,drop=F]
  v[v<0]=0
  return(cbind(res,v))
}


##########################################################
.GRmsPeakSimple3<-function (y, noise.local = NULL, span = 5, snr.thresh = 1,compCl=T) {
  
  nvar <- length(y)
  index <- .GRmsExtrema(y, span = span)
  
  index.min = index$index.min
  index.max = index$index.max
  
  pksnull = data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("tick.loc","tick.left","tick.right","SNR","aboveSN"))), stringsAsFactors=F)
  
  imax <- which(index.max)
  tick.loc <- imax#[good.snr]
  
  if ((npeak = length(tick.loc)) == 0) return(pksnull)
  tick.left <- tick.right <- rep(1, length(tick.loc))
  for (i in 1:length(tick.loc)) {
    #for (j in tick.loc[i]:1) 
    tick.left[i] <-max(which(index.min[1:tick.loc[i]]))
    #for (j in tick.loc[i]:length(x)) 
    tick.right[i] <- min(which(index.min[tick.loc[i]:nvar]))+tick.loc[i]-1
  }
  if (tick.right[length(tick.loc)] == 1) tick.right[length(tick.loc)] <- nvar
  if (tick.left[1] == 1) tick.left[1] <- which.min(y[1:tick.loc[1]]!=0)
  
  #    keep <- !duplicated(tick.left)
  #    tick.loc <- tick.loc[keep]
  #    tick.left <- tick.left[keep]
  #    tick.right <- tick.right[keep]
  for(ix in 1:length(tick.loc)){
    l=tick.left[ix]:tick.right[ix]
    tick.left[ix]=max(tick.left[ix],min(l[which(y[l]>0)])-1)
    tick.right[ix]=min(tick.right[ix],max(l[which(y[l]>0)])+1)
  }
  ###  check duplicated left -> use the highest peak
  lpks=lapply(unique(tick.left),function(x) which(tick.left==x))
  pks=do.call("rbind",lapply(lpks,function(ix){
    if(length(ix)>1) return(c(tick.loc=tick.loc[ix[which.max(y[tick.loc[ix]])]],tick.left=min(tick.left[ix]), tick.right=max(tick.right[ix])))
    c(tick.loc=tick.loc[ix],tick.left=tick.left[ix], tick.right=tick.right[ix])
  }))
  
  
  if(!is.null(noise.local)) snr <- y[imax]/noise.local[imax] else snr <- y[imax]
  good.snr <- (snr > snr.thresh)
  snr <- snr[good.snr]
  abosn=imax[good.snr]
  pks=cbind(pks,SNR=y[pks[,1]],aboveSN=pks[,1]%in%abosn)
  
  if(compCl){
    icl=0;lcl=l=c()
    for(ix in 1:nrow(pks)){
      if(any((pks[ix,3]:pks[ix,2])%in%l)){
        l=c(l,pks[ix,3]:pks[ix,2])
      }
      else{
        l=pks[ix,3]:pks[ix,2]
        icl=icl+1      
      }
      lcl=c(lcl,icl)    
    }
    pks=cbind(pks,cl=lcl)
  }
  
  return(pks)
  
}

.GRgroupPks<-function(apks,tabres){
 
  inloc=inlocr=inlocl=rep(0,nrow(apks))
  lseg=sort(unique(tabres[tabres[,4]!=0,2]))
  lseg=lseg[lseg>0]
  lsegsc=lapply(lseg,function(x) which(tabres[,2]==x))
  for(ij in lseg) inloc[apks$tick.loc%in%which(tabres[,2]==ij)]=ij
  for(ij in lseg) inlocl[apks$tick.left%in%which(tabres[,2]==ij)]=ij
  for(ij in lseg) inlocr[apks$tick.right%in%which(tabres[,2]==ij)]=ij
  toadd=cbind(pkreg=inloc,lefreg=inlocl,rireg=inlocr)
  ### remove all zeros alone in their cluster
  allz=apply(toadd,1,function(x) all(x==0))*1
  scl=paste(apks$samp,apks$cl)
  tab=table(scl,factor(allz,levels=0:1))
  l2rm=which(scl%in%names(which(tab[,1]==0 & tab[,2]>0)))
  if(length(l2rm)){toadd=toadd[-l2rm,,drop=F];apks=apks[-l2rm,,drop=F]}
  if(nrow(apks)==0) return(NULL)
  ### chk if all agree
  ndif=apply(toadd,1,function(x) length(unique(x[x>0])))
  ndif[which(ndif>1)]=NA
  ndif[which(ndif==1)]=apply(toadd[which(ndif==1),,drop=F],1,function(x) unique(x[x>0]))
  if(length(ndif)>1) if(any(ndif==0,na.rm=T)) for(mergezer in 1:2){
    ### all zeros before
    if(any(ndif[-length(ndif)]==0,na.rm=T)){
      lzer=which(ndif[-length(ndif)]==0)
      lzer=lzer[which(apks$samp[lzer+1]==apks$samp[lzer] & apks$cl[lzer+1]==apks$cl[lzer] & ndif[lzer+1]>0)]
      if(length(lzer)>0) ndif[lzer]=ndif[lzer+1]
    }
    ### all zeros after
   if(any(ndif[-1]==0,na.rm=T)){
      lzer=which(ndif[-1]==0)+1
      lzer=lzer[which(apks$samp[lzer-1]==apks$samp[lzer] & apks$cl[lzer-1]==apks$cl[lzer] & ndif[lzer-1]>0)]
      if(length(lzer)>0) ndif[lzer]=ndif[lzer-1]
    }
  }
  ### chk if only in the sample/cluster NA
  if(any(is.na(ndif))){
    scl=paste(apks$samp,apks$cl)
    lna=which(is.na(ndif))
    l=which(scl%in%scl[lna])
    if(!all(is.na(ndif[l]))){
    tab=table(scl[l],ndif[l])
    tab=tab[rowSums(tab>0)==1,,drop=F]
    if(nrow(tab)>0) for(i in rownames(tab)) ndif[which(scl%in%i)]=unique(na.omit(ndif[which(scl%in%i)]))
    }
  }  
  
  ### check maximum overlap
  if(any(is.na(ndif))){
    for(i in which(is.na(ndif))){
    llpk=apks$tick.left[i]:apks$tick.right[i]
    maxov=which.max(sapply(lsegsc,function(x) length(intersect(x,llpk))/length(llpk)))
    if(maxov%in%toadd[i,]) ndif[i]=maxov
    }
  }
  return(cbind(apks,toadd,Pk=ndif))
}


#####################################################################
###
.oGRgroupPks<-function(apks,lpks=NULL,lseg=rep(1,length(lpks))){
  
  
  apks0=apks
  apks=apks[(apks$pkreg+apks$lefreg+apks$rireg)>0,,drop=F]
  if(nrow(apks)==0) return(NULL)
  regs=apks[,c("pkreg","lefreg","rireg"),drop=F]
  lregs=unique(as.vector(as.matrix(regs)))
  lregs=lregs[lregs>0]
  uregs=as.character(regs[,"pkreg"])
  l1=which(uregs==0 & (regs[,"rireg"]==0 | regs[,"lefreg"]==regs[,"rireg"]))
  if(length(l1)>0) uregs[l1]=regs[l1,"lefreg"]
  l1=which(uregs==0 & regs[,"lefreg"]==0)
  if(length(l1)>0) uregs[l1]=regs[l1,"rireg"]
  
  l1=which(uregs==0 & (regs[,"lefreg"]!=regs[,"rireg"]))
  uregs[l1]=NA
  ############
  # Merging 
  uuregs=paste(uregs,"1",sep=".")
  if(any(table(lseg)>1)){
    for(j in names(which(table(lseg)>1))){
      l=which(uregs==j)
      ipks=lpks[which(lseg==j)]
      uuregs[l]=unlist(tapply(l,apks$samp[l],function(lx) names(ipks[.GRchk2vec(apks$tick.loc[lx],ipks)])))
    }
  }
  uregss=paste(apks$samp,uuregs,sep=".")
  apks$uuregs=uuregs
  ###

  ldups=names(which(table(uregss)>1))
  napks=apks
  if(length(ldups)>0){
    napks=list()
    for(i in ldups){
      l=which(uregss==i)
      lmax=l[which.max(apks$SNR[l])]
      t2add=apks[lmax,]
      t2add$tick.left=min(apks$tick.left[l])
      t2add$tick.right=max(apks$tick.right[l])
      napks[[i]]=t2add
    }
    napks=rbind(apks[which(!uregss%in%ldups),,drop=F],do.call("rbind",napks))
    uregss=c(uregss[which(!uregss%in%ldups)],ldups)
  }
  napks=napks[order(napks$samp,napks$tick.loc),,drop=F]
  
  ###########
  #check if NAs
#   lnas=grep("NA",rownames(napks))
#   if(length(lnas)>0){
#     for(ina in lnas){
#       isamp=napks[ina,1]
#       isc=napks[ina,2]
#       li=which(napks[,1]==isamp);
#       
#       liu=li[napks[li,4]<napks[ina,2]]
#       liu=liu[which.min(abs(napks[liu,4]-isc))]
#       napks[liu,4]=napks[ina,2]
#       liu=li[napks[li,3]>napks[ina,2]]
#       liu=liu[which.min(abs(napks[liu,3]-isc))]
#       napks[liu,3]=napks[ina,2]
#     }
#     napks=napks[-lnas,,drop=F]
#   }
  return(napks)
}


###
# merge 2 vec y on x
.GRchk2vec<-function(x,y){
  nx=length(x)
  ny=length(y)
  if(nx==1) return(which.min(abs(y-x)))
  if(ny==1){lma=rep(NA,length(x));lma[which.min(abs(x-y))]=1;return(lma)}
  if(nx==ny) return(order(x))
  if(nx>(ny+5)){
    lma=rep(NA,nx)
    for(ic in 1:ny){
      imi=which(is.na(lma))[which.min((y[ic]-x[which(is.na(lma))])^2)]
      lma[imi]=ic
    }
    if(any(is.na(lma)))
      for(ic in which(is.na(lma))) lma[ic]=which.min((y-x[ic])^2)
      
      return(lma)
  }
  
  if(ny<nx){ ### more pekas than in the model peak
    mc=(combn(nx,ny)) 
    lma=mc[,which.min(apply(mc,2,function(ix) sum((y-x[ix]))^2))]
    lma=match(1:nx,lma)
    l1=which(is.na(lma))
    l1=l1[l1<which(!is.na(lma))[1]]
    if(length(l1)>0) lma[l1]=lma[which(!is.na(lma))[1]]
    l1=which(is.na(lma))
    l1=l1[l1>rev(which(!is.na(lma)))[1]]
    if(length(l1)>0) lma[l1]=lma[rev(which(!is.na(lma)))[1]]
    if(any(is.na(lma)))
      for(ic in which(is.na(lma))) lma[ic]=which.min((y-x[ic])^2)
    
    return(lma)
  }
  if(ny>nx){
    mc=combn(ny,nx)
    return(mc[,which.min(apply(mc,2,function(ix) sum((x-ix))^2))])
  }
}
