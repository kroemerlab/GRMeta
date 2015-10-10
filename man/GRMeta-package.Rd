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

\dontrun{
library(GRMeta)
library(chron)

(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
lcdata=loadAgilentData(ifile =testfiles[1],params = list(AssayName="lc"))
gcdata1=loadAgilentData(ifile =testfiles[2],params = list(AssayName="gc"))
gcdata2=loadAgilentData(ifile =testfiles[3],params = list(AssayName="gc"))

finaldata=mergeSet(lcdata,gcdata1,gcdata2)

## or step by step
gcdata=mergeSet(gcdata1,gcdata2)
gcdata=remove(gcdata,Sid=c("Sa_cond1_5.B2"))
gcdata=update(gcdata,what="Sid",formerid = c("Sa_cond1_5.B1"),newid = c("Sa_cond1_5"))

finaldata=mergeSet(lcdata,gcdata)



### Export "Area", replace NA by "-" and "NA"
exportDataXLSX(mydata,outfile="test.xlsx",ldata="Area",sortSid="Sid",characterNA=c("-","NA"))
}


}
