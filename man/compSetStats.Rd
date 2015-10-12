\name{compSetStats}
\alias{compSetStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute sample statistics}
\description{
Compute sample statistics in a metaboSet object}
\usage{
compSetStats(obj, lStats = list(), whatSamp = "Sa", whatQC = "QC")
}
\arguments{
  \item{obj}{metaboSet object}
  \item{lStats}{List of stats to be performed}
  \item{whatSamp}{sType level or vector of Sids for samples}
  \item{whatQC}{sType level or vector of Sids for QCs}
}
\details{
}
\value{Matrix of stats}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
mydata=loadAgilentData(ifile =testfiles[1])

##### Number of "Sa" with Height>0,  Height>100 | Proportion with Height>0,  Height>100  | Proportion of non-NA with Height>0,  Height>100 
lStats=list(c("Sum","Sa","Height>0"),c("Sum","Sa","Height>100"),c("Prop","Sa","Height>0"),c("Prop","Sa","Height>100"),
c("PropNNA","Sa","Height>0"),c("PropNNA","Sa","Height>100"))
result=compSetStats(mydata,lStats=lStats,whatSamp="Sa")
head(result)

##### Same results if actual Sid are given
lSamp=mydata$Sid[mydata$Meta$sType=="Sa"]
result=compSetStats(mydata,lStats=lStats,whatSamp=lSamp)
head(result)


##### Compute CVs on QC and Sa with at least 1 non NAs or at least 5
result=compSetStats(mydata,lStats=list(c("CV","QC","Area",1),c("CV","Sa","Area",1),c("CV","QC","Area",5),c("CV","Sa","Area",5)))
head(result)



}
\keyword{misc}
