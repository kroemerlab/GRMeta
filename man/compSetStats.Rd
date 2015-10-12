\name{compSetStats}
\alias{compSetStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute sample statistics}
\description{
Compute sample statistics in a metaboSet object}
\usage{
compSetStats(obj, lStats = list(), whatSamp = "Sa", whatQC = "QC",add=FALSE)
}
\arguments{
  \item{obj}{metaboSet object}
  \item{lStats}{List of stats to be performed}
  \item{whatSamp}{sType level or vector of Sids for samples}
  \item{whatQC}{sType level or vector of Sids for QCs}
  \item{add}{If TRUE the stats will be added to the obj$Annot and the result object returned}
}
\details{
}
\value{Matrix of stats or metaboSet object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
mydata=loadAgilentData(ifile =testfiles[1])
mydata

##### Here "St" will be treated as sample
##### Number of "St" with Height>0,  Height>100 | Proportion with Height>0,  Height>100  | Proportion of non-NA with Height>0,  Height>100 
lStats=list(c("Sum","Sa","Height>0"),c("Sum","Sa","Height>100"),c("Prop","Sa","Height>0"),c("Prop","Sa","Height>100"),
c("PropNNA","Sa","Height>0"),c("PropNNA","Sa","Height>100"))
result=compSetStats(mydata,lStats=lStats,whatSamp="St")
head(result)

##### Same results if actual Sid are given
lSamp=mydata$Sid[mydata$Meta$sType=="St"]
result=compSetStats(mydata,lStats=lStats,whatSamp=lSamp)
head(result)


##### Compute CVs on QC and Sa with at least 1 non NAs or at least 5
result=compSetStats(mydata,lStats=list(c("CV","QC","Area",1),c("CV","Sa","Area",1),c("CV","QC","Area",5),c("CV","Sa","Area",5)))
head(result)

##### Compute CVs on QC with at least 2 non NAs / sum of SBR>5 on "St" / add to the object
mydata=addSBRMatrix(mydata,thresh=10,SidBlank = c("Bl_prerun-eau_1","Bl_prerun-eau_1b"),newname = "SBR")
result=compSetStats(mydata,lStats=list(c("CV","QC","Area",2),c("Sum","Sa","SBR>5")),whatSamp="St",add=TRUE)
head(result$Annot)




}
\keyword{misc}
