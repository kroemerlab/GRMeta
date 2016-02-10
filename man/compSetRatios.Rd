\name{compSetRatios}
\alias{compSetRatios}
\title{
Compute metabolite ratios}
\description{
Compute metabolite ratios from a metaboSet object}
\usage{
compSetRatios(obj, l2add = names(obj$Data), lpairs = RatiosGR, conv2set = FALSE)
}
\arguments{
  \item{obj}{metaboSet object}
  \item{l2add}{List of matrices in obj for which the ratios are computed}
  \item{lpairs}{List of pairs: if it loads \code{\link{RatiosGR}} supported}
  \item{conv2set}{Should the resulting calculation converted to metaboSet object}
}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
mydata1=loadAgilentData(ifile =testfiles[1],params = list(AssayName="assay1"))
res=compSetRatios(mydata1,l2add="Area",conv2set=FALSE)
res

## convert to metaboSet
res=compSetRatios(mydata1,l2add="Area",conv2set=TRUE)
res

## 
mydata2=loadAgilentData(ifile =testfiles[2],params = list(AssayName="assay2"))
cmbdata=mergeSet(mydata1,mydata2)
res=compSetRatios(cmbdata,l2add="Area",conv2set=TRUE)
res

}
\keyword{ misc }
