\name{GRMeta-package}
\alias{GRMeta-package}
\alias{GRMeta}
\docType{package}
\title{
GRMeta
}
\description{
GRMeta is great
}
\details{
The DESCRIPTION file:
}
\author{
GRMeta
}
\keyword{ package }
\examples{

library(GRMeta)

(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
lcdata=loadAgilentData(ifile =testfiles[1],params = list(AssayName="lc"))
gcdata1=loadAgilentData(ifile =testfiles[2],params = list(AssayName="gc",Batch="B1"))
gcdata2=loadAgilentData(ifile =testfiles[3],params = list(AssayName="gc",Batch="B2"))

finaldata=mergeSet(lcdata,gcdata1,gcdata2)
finaldata
## or step by step
gcdata=mergeSet(gcdata1,gcdata2)
gcdata=deleteEntries(gcdata,Sid=c("Sa_cond1_5.B2"))
gcdata=update(gcdata,what="Sid",formerid = c("Sa_cond1_5.B1"),newid = c("Sa_cond1_5"))

finaldata=mergeSet(lcdata,gcdata)

## plot analytes that correspond to AMP
plot(lcdata,MetName="AMP")

## plot all analytes from method="lc"
plot(finaldata,Method="lc")

## by default all analytes are considered if dumping to PDF
\dontrun{
plot(finaldata,outfile="All.pdf")
}

### Export "Area", replace NA by "-" and "NA"
\dontrun{
exportSetXLSX(mydata,outfile="test.xlsx",ldata="Area",sortSid="Sid",characterNA=c("-","NA"))
}


}
