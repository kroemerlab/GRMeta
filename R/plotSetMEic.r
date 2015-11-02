
plotSetMEic<-function(obj,WhichRT="rtcor",groupCol=NULL,colorCol=NULL,addfile="./",endfile="-Eic",repDots="-",cexEL=0.6,...){

  
  dots<-list(...)
  
  if(is.null(obj$Eic)) stop("Object does not contain any EIC infos.")
  
  lanalytes=NULL
  for(k in names(obj$Eic$File))  if(!is.null(dots[[k]])){
    l2add=obj$Eic$File$Analyte[obj$Eic$File[,k]%in%dots[[k]]]
    l2add=l2add[!l2add%in%lanalytes]
    if(length(l2add)>0) lanalytes=c(lanalytes,l2add)
  }
  if(is.null(lanalytes)){
    lanalytes=obj$Eic$File$Analyte
  }
  
  lanalytes=lanalytes[lanalytes%in%obj$Eic$File$Analyte]
  lfiles=lfiles2=unique(obj$Eic$File$EicFile[obj$Eic$File$Analyte%in%lanalytes])
  if(!is.null(obj$Eic$Path)) lfiles2=paste(obj$Eic$Path,lfiles2,sep="")
  lexist=which(file.exists(lfiles2))
  lfiles=lfiles[lexist]
  lanalytes=lanalytes[lanalytes%in%obj$Eic$File$Analyte[obj$Eic$File$EicFile%in%lfiles]]
  lfiles=lfiles[lfiles%in%obj$Eic$File$EicFile]
  cat(length(lanalytes),"analytes from",length(lfiles),"files\n")
  if(length(lanalytes)==0) stop("Nothing found!")
  
  width = ifelse(is.null(dots$width),5*par()$mfrow[2],dots$width)
  height = ifelse(is.null(dots$height),5*par()$mfrow[1],dots$height)
  dots=dots[names(dots)%in%names(par())]
  
  
  grp=obj$Meta$sType
  if(!is.null(groupCol)) if(!is.na(groupCol)) if(colorCol%in%names(obj$Meta)) grp=obj$Meta[,groupCol]
  grp=factor(grp)
  llsids=tapply(obj$Sid,grp,c)
  
  cols=brewer.pal(9,"Set1")[-6][as.numeric(factor(obj$Meta$sType))]
  if(!is.null(colorCol)) if(!is.na(colorCol)) if(colorCol%in%names(obj$Meta)) cols=obj$Meta[,colorCol]
  cols[is.na(cols)]="black"
  names(cols)=obj$Sid
  
  par.def=par(no.readonly = TRUE)
  
for(ifeic in lfiles){

l2plots=obj$Eic$File$Analyte[obj$Eic$File$EicFile==ifeic]
l2plots=l2plots[l2plots%in%lanalytes]
if(!is.null(obj$Eic$Path)) ifeic=paste(obj$Eic$Path,ifeic,sep="")
cat("Found",ifeic)
load(ifeic)
cat(": ")

lx=round(seq(1,length(l2plots),length.out = 5)[2:4])
# if(sepx[2]%in%names(idf)) llsids=tapply(idf$Sid,idf[,sepx[2]],unique)
#for(ix in 1:length(l2plots)){
  for(ix in 1:length(l2plots)){
    cat(ifelse(ix%in%lx,"X","."))
    analyte=l2plots[ix]
ieicpk=obj$Eic$File[analyte,]$EicPK
ieic=obj$Eic$File[analyte,]$EicId

whichrt=WhichRT
if(is.null(whichrt))  whichrt="rtcor"
if(!whichrt%in%names(dfeic)) whichrt="rt"

ceic=dfeic[dfeic$Eic==ieic & dfeic$Samp%in%obj$Eic$Sample$Samp,]
ceic$Sid=obj$Eic$Sample$Sid[match(ceic$Samp,obj$Eic$Sample$Samp)]
ceic$cols=cols[ceic$Sid]
rtr=range(ceic[,whichrt])
rtr=rtr+c(-.5,1)*0.1*diff(rtr)
rtr[1]=max(0,rtr[1])
rtr=range(pretty(seq(rtr[1],rtr[2],length.out = 7)))


ipkmat=eicpk[[ieic]]$Pks[[ieicpk]]
ipkmat=ipkmat[rownames(ipkmat)%in%obj$Eic$Sample$Samp,]
rownames(ipkmat)=obj$Eic$Sample$Sid[match(rownames(ipkmat),obj$Eic$Sample$Samp)]
Mint=obj$Eic$File[analyte,]$Mint

outfile=paste(addfile,gsub("\\.",repDots,analyte),endfile,".pdf",sep="")
pdf(file=outfile,width=width,height=height)
par(dots)
.plotEIC(ceic,ipkmat,whichrt,llsids,Mint=Mint,rtr=rtr,cexEL=cexEL)
dev.off()

}
cat("\n")
} ## end of for ifeic
  on.exit(par(par.def))
}
