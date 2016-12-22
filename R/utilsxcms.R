## need xcms installed


######## rewrite of findROI (https://github.com/nathaniel-mahieu/centWaveP)

.GRfindROIs <- function(xr, ppm = 2, prefilter = c(0,0), maxskip = 0) {
  #SEXP findmzROI(SEXP mz, SEXP intensity, SEXP scanindex, SEXP mzrange, SEXP scanrange, SEXP lastscan, SEXP dev, SEXP minEntries, SEXP prefilter, SEXP noise, SEXP maxskip)
  .Call("findmzROI",
        as.double(xr@env$mz),
        as.double(xr@env$intensity),
        as.integer(xr@scanindex), 
        as.double(range(xr@env$mz)),
        as.integer(c(1,length(xr@scantime))), 
        as.integer(length(xr@scantime)),
        as.double(ppm * 1E-6), 
        as.integer(0), 
        as.integer(prefilter), 
        as.integer(0), 
        as.integer(maxskip),
        PACKAGE="GRMeta")
}

######### expanded xcms:::rawMat

.GRrawMat<-function (object, mzrange = numeric(), rtrange = numeric(), 
                   scanrange = numeric(),convsec=T,padsc=FALSE,naVal=NA) {
  require(xcms)
  
  if (length(rtrange) >= 2) {
    rtrange <- range(rtrange)
    scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= 
                                                    rtrange[2])
    scanrange <- c(match(TRUE, scanidx), length(scanidx) - 
                     match(TRUE, rev(scanidx)))
  }
  else if (length(scanrange) < 2) 
    scanrange <- c(1, length(object@scantime))
  else scanrange <- range(scanrange)
  startidx <- object@scanindex[scanrange[1]] + 1
  endidx <- length(object@env$mz)
  if (scanrange[2] < length(object@scanindex)) 
    endidx <- object@scanindex[scanrange[2] + 1]
  #  print(scanrange)
  scans <- rep(scanrange[1]:scanrange[2], diff(c(object@scanindex[scanrange[1]:scanrange[2]], endidx)))
  rtrange <- c(object@scantime[scanrange[1]], object@scantime[scanrange[2]])
  masses <- object@env$mz[startidx:endidx]
  ids=startidx:endidx
  int <- object@env$intensity[startidx:endidx]
  massidx <- 1:length(masses)
  if (length(mzrange) >= 2) {
    mzrange <- range(mzrange)
    massidx <- massidx[(masses >= mzrange[1]) & (masses <=  mzrange[2])]
  }
  else mzrange <- range(masses)
  y <- int[massidx]
  ret=cbind(scan=scans[massidx],mz = masses[massidx],y=y, id=ids[massidx])
  if(padsc){
    lusc=scanrange[1]:scanrange[2]
    lusc=lusc[!lusc%in%ret[,"scan"]]
    if(length(lusc)) ret=rbind(ret,cbind(scan=lusc,mz=NA,y=NA,id=NA))
  }
  ret=cbind(ret,rt = object@scantime[ret[,"scan"]])
  if(convsec) ret[,"rt"]=ret[,"rt"]/60
  if(!is.na(naVal)) ret[which(is.na(ret[,"y"]) | ret[,"y"]<naVal),"y"]=naVal
  ret[order(ret[,"scan"],ret[,"mz"]),,drop=F]
}

########## compute SBR from Agilent cefmat data.frame
.GRcompSBRCef<-function(blsamp,cefmat,dppm=21,drt=.5,nSmo=31,minNoise=350){
  
  require(xcms)
  xRawBl <- xcmsRaw(blsamp, profmethod = "bin", profparam = list(step = 0.001), profstep = 0)
  rangeSc=c(1,length(xRawBl@scantime))
  dscan=ceiling(60*drt/median(diff(xRawBl@scantime)))+nSmo
  
  strt0=Sys.time()
  aybs=rep(minNoise,nrow(cefmat))
  for(ientry in 1:nrow(cefmat)){
    x=cefmat[ientry,]
    lmz=range(x$mz*(1+c(-1,1)*dppm*10^-6))
    iscan=which.min(abs(xRawBl@scantime-x$rt*60))
    lsc=.GRchkrange(iscan+c(-1,1)*dscan,rangeSc)
    
    xb=xcms:::rawEIC(xRawBl, mzrange =lmz, scanrange = lsc)
    if(!any(xb$intensity<minNoise)) next
    xbi=xb$intensity[match(lsc[1]:lsc[2],xb$scan)]
    xbi[xbi<minNoise]=minNoise
    xbi=.GRdoksmooth(list(x=lsc[1]:lsc[2],y=xbi),nSmo)
    xb$ybs=round(approx(xbi$x,xbi$y,xb$scan)$y,1)
    aybs[ientry]=xb$ybs[xb$scan==iscan]
  }
  names(aybs)=cefmat$mfid2
  invisible(aybs)
}
