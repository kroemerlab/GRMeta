
xcms2mzXML<-function(xr,filename,lscan2export=NULL,comment="GRMeta"){

  if("character"%in%class(xr)){
    if(!file.exists(xr)) stop(paste0(xr," does not exists!"))
    xr=xcmsRaw(xr)
  }
  infile=unclass(xr@filepath)
  
  if(is.null(lscan2export)) lscan2export=1:length(xr@scantime)
  lscan2export=lscan2export[lscan2export%in%(1:length(xr@scantime))]
  rtrange=range(xr@scantime[lscan2export])
  
  
fprintf = function(fp, level, ..., append=TRUE)
{ # helper function
  x = paste(..., sep="")
  if (length(x)==0 || is.null(x)) return(NULL)
  spaces = if (level>0) paste0(rep("  ", level)) else ""
  x = gsub("'", "\"", x)
  cat(spaces, x, file=fp, sep="")
  NULL
}  # done with local functions

#### Form header
mzXML=list()
Str = "http://sashimi.sourceforge.net/schema_revision/mzXML_3.1"
mzXML$header = paste0( "<mzXML xmlns='",Str,"'\n  ",
                       "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n  ",
                       "xsi:schemaLocation='", Str, " ", Str, "/mzXML_idx_3.1.xsd'>\n")
mzXML$parentFile =paste0( "  <parentFile filename='",infile,"' ",
                          "fileType='RAWData' fileSha1='0000000000000000000000000000000000000000'/>\n")
Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
mzXML$dataProcessing = paste0("  <dataProcessing>\n",
                              "   <software type='processing' name='",comment,"' version='1.0.0.0' completionTime='",Time,
                              "'/>\n  </dataProcessing>\n")
if(file.exists(infile)){
xr2=mzR:::openMSfile(infile)
msInfos=instrumentInfo(xr2)
rm(list='xr2')
} else msInfos=list(manufacturer="unknown",model="unknown",ionisation="unknown",analyzer="unknown")


mzXML$msInstrument = paste0("  <msInstrument>\n",
                            "   <msManufacturer category='msManufacturer' value='",msInfos$manufacturer,"' />\n",
                            "   <msModel category='msModel' value='",msInfos$model,"' />\n",
                            "   <msIonisation category='msIonisation' value='",msInfos$ionisation,"' />\n",
                            "   <msMassAnalyzer category='msMassAnalyzer' value='",msInfos$analyzer,"' />\n",
                            "   <software type='acquisition' name='Xcalibur' version='unknown' />\n","  </msInstrument>\n")


cat("Exporting ",length(lscan2export)," scans to '",filename,"' ",sep="")

fp  = file(filename, "w")
fprintf(fp, 0, "<?xml version='1.0' encoding='ISO-8859-1'?>\n", append=FALSE)
fprintf(fp, 0, mzXML$header)
fprintf(fp, 0, sprintf(" <msRun scanCount='%d' startTime='PT%.5fS' endTime='PT%.5fS'>\n",length(lscan2export),rtrange[1],rtrange[2]))
fprintf(fp, 0, mzXML$parentFile)
fprintf(fp, 0, mzXML$msInstrument)
fprintf(fp, 0, mzXML$dataProcessing)


ix=1
indexScan = list()
lperc=round(seq(1,length(lscan2export),length.out = 10))
for(idx in 1:length(lscan2export)){
if(idx%in%lperc) cat("x")
  indexScan[[idx]]=paste0("  <offset id='",idx,"'>",seek(fp),"</offset>\n")
  ix=lscan2export[idx]
  v=c( "msLevel"=1)
  p=xcms:::getScan(xr,ix)
  v['peaksCount']=nrow(p)
  v['polarity']=ifelse(xr@polarity[ix]=="negative","-","+")
  v['scanType']='Full'
  v['retentionTime']=sprintf("PT%.5fS",xr@scantime[ix])
  v['lowMz']=sprintf("%.5f",min(p[,1]))
  v['highMz']=sprintf("%.5f",max(p[,1]))
  bp=p[which.max(p[,2]),]
  v['basePeakMz']=sprintf("%.5f",bp[1])
  v['basePeakIntensity']=sprintf("%e",round(bp[2]))
  v['totIonCurrent']=sprintf("%e",round(sum(p[,2])))
  
  
ScanHeader=paste0("   ",names(v),"='",v,"'")
ScanHeader[length(ScanHeader)]=paste0(ScanHeader[length(ScanHeader)],">\n")
ScanHeader=paste(c(paste0("  <scan num='",idx,"'"),ScanHeader),collapse="\n")

p=paste0("   <peaks precision='32'\n    byteOrder='network'\n    contentType='m/z-int'\n    compressionType='none'\n    compressedLen='0' >",
         caTools::base64encode(as.vector(t(p)), endian="big", size=4), "</peaks>\n")
fprintf(fp, 0, ScanHeader)
fprintf(fp, 0, p)
fprintf(fp, 0, "  </scan>\n")
}

fprintf(fp, 0, " </msRun>\n")
cat("\n")

#### write offset
n = seek(fp)
fprintf(fp, 0, " <index name='scan'>\n")
fprintf(fp, 0, paste(indexScan,collapse=""))
fprintf(fp, 0, " </index>\n")
fprintf(fp, 0, " <indexOffset>",n,"</indexOffset>\n")
cat(" <sha1>", file=fp, sep="")
n = seek(fp)
close(fp)
sha1 = digest::digest(filename, algo="sha1", file=TRUE, length=n)
cat(sha1, "</sha1>\n</mzXML>\n", file=filename, append=TRUE, sep="")
invisible(NULL)
}
