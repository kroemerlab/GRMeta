############################
.GRBlankEIC1<-function(xRawBl=NULL,lmz,lsc,minNoise,bsllambda){
  
  lsc=min(lsc):max(lsc)
  if(is.null(xRawBl)) return(rep(minNoise,length(lsc)))
#  lscbl=.GRchkrange(lsc,rangeScBl[c(1,4)])
  xb=xcms:::rawEIC(xRawBl, mzrange =lmz)
  xbc=xb$scan
  xbi=xb$intensity
  xbi=xbi[match(lsc,xbc)]
  xbi[xbi<minNoise | is.na(xbi)]=minNoise
  #   xbi=.GRdoksmooth(list(x=lscbl[1]:lscbl[2],y=xbi),parDeco$nSmo2)
  #   eic[,"ybl"]=round(approx(xbi$x,xbi$y,eic[,"scan"])$y,1)*parDeco$sbr
  ybl=.GRasysm(xbi,p=0.5,lambda = bsllambda)  ## use 50% as noise baseline
  ybl[ybl<=minNoise | is.na(ybl)]=minNoise
  ybl
}
############################
.GRingetROIBl<-function(xRawBl,lmz,lsc,rangeScBl,minNoise,bsllambda){
  
  lscbl0=min(lsc):max(lsc)
  lscbl=.GRchkrange(lsc,rangeScBl[c(1,4)])
  xb=xcms:::rawEIC(xRawBl, mzrange =lmz, scanrange = lscbl)
  xbc=xb$scan
  xbi=xb$intensity
  lscbl=range(c(lscbl,lsc))
  xbi=xbi[match(lscbl[1]:lscbl[2],xbc)]
  xbi[xbi<minNoise | is.na(xbi)]=minNoise
  #   xbi=.GRdoksmooth(list(x=lscbl[1]:lscbl[2],y=xbi),parDeco$nSmo2)
  #   eic[,"ybl"]=round(approx(xbi$x,xbi$y,eic[,"scan"])$y,1)*parDeco$sbr
  ybl=.GRasysm(xbi,p=0.5,lambda = bsllambda)  ## use 50% as noise baseline
  ybl=ybl[match(lscbl0,min(xbc):max(xbc))]
  ybl[ybl<=minNoise | is.na(ybl)]=minNoise
  ybl
}

##############################
.GRgetEICroi<-function(xroi,xRaw,parDeco,addSc=floor(parDeco$addSc/2)){
  
  lmz=xroi[2:3];lsc=xroi[4:5]
  lsc2=.GRchkrange(lsc+c(-1,1)*addSc,parDeco$rangeSc[c(1,4)])
  
  xs=xcms:::rawEIC(xRaw, mzrange = lmz, scanrange =lsc2)
  xs$mz=xcms:::rawMZ(xRaw,mzrange = lmz, scanrange = lsc2)
  xs$mz[xs$mz==0]=NA;xs$intensity[xs$intensity<parDeco$minZero]=parDeco$minZero
  xs$y=xs$intensity
  xs$inwin=(xs$scan%in%lsc[1]:lsc[2])*1
  eic=do.call("cbind",xs[names(xs)!="intensity"])
 
  ##### smooth -> include gap filling here
  ys=signal::sgolayfilt(xs$y, p = 2, n = parDeco$sg.window)
  ys[ys<parDeco$minZero]=parDeco$minZero
  eic=cbind(eic,ys=ys)
  
  ### compute noise
  tmpeic=eic[,"y"] %>% {.[.==0] = eic[.==0,"ys"]*0.5; .}
  rollav=stats:::filter(eic[,"y"], 1/rep(parDeco$sg.window, parDeco$sg.window)) %>% { .[is.na(.)] = eic[,"y"][is.na(.)]; . }
  noise.sd2 <- caTools::runsd(tmpeic-rollav,parDeco$sg.window)
  eic=cbind(eic,ys.sd=noise.sd2)
  
  ### compute baseline
  re=.GRbslrf(eic[,"scan"],eic[,"ys"],NoXP = parDeco$addSc)
  bsl<-re$fit
  bsl[bsl<=parDeco$minZero]=parDeco$minZero
  
  ys.score<- (eic[,"ys"]-bsl)/max(re$sigma,10^-3)
  ys.score[which(abs(ys.score)>10^5)]=sign(ys.score[which(abs(ys.score)>10^5)])*10^5
  #bsl<-.GRasysm(eic[,"y"],p=bsl.p,lambda = parDeco$bsl.lambda)
  eic=cbind(eic,bsl=bsl,ys.score=ys.score)
  return(eic)
}

