\name{compSetStats}
\alias{compSetStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute sample statistics}
\description{
Compute sample statistics in a metaboSet object}
\usage{
compSetStats(obj, lStats = list(),sType=NULL,Sids=NULL,addtoname=NULL,addtoobject=FALSE)
}
\arguments{
  \item{obj}{metaboSet object}
  \item{lStats}{List of stats to be performed}
  \item{sType}{sType level for samples}
  \item{Sids}{List Sids for the calculation}
  \item{addtoname}{Concatenate string to the begining of the stat names}
  \item{addtoobject}{If TRUE the stats will be added to the obj$Annot and the result object returned}
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

##### Use sample of sType=="St"
##### Number with Height>0,  Height>100 | Proportion with Height>0,  Height>100  | Proportion of non-NA with Height>0,  Height>100 
lStats=list(c("Sum","Height>0"),c("Sum","Height>100"),c("Prop","Height>0"),c("Prop","Height>100"),
c("PropNNA","Height>0"),c("PropNNA","Height>100"))
result=compSetStats(mydata,lStats=lStats,sType="St")
head(result)

##### Same results if actual Sid are given
lSamp=mydata$Sid[mydata$Meta$sType=="St"]
result=compSetStats(mydata,lStats=lStats,Sids=lSamp)
head(result)


##### Compute CVs on QCwith at least 1 non NAs or at least 5 and add "myQC" before the column name
result=compSetStats(mydata,lStats=list(c("CV","Area",1),c("CV","Area",5)),sType="QC",addtoname="myQC")
head(result)

##### Compute CVs on QC with at least 2 non NAs and check the number of "St" with SBR>5 / add to the object
mydata2=addSBRMatrix(mydata,thresh=10,SidBlank = c("Bl_prerun-eau_1","Bl_prerun-eau_1b"),newname = "SBR")
## Compute CVs on QC with at least 2 non NAs
mydata2=compSetStats(mydata2,lStats=list(c("CV","Area",2)),sType="QC",addtoname="myQC",addtoobject=TRUE)
## Compute on St: CV with at least 2 non NAs, number of SBR>5 and number of non missing
mydata2=compSetStats(mydata2,lStats=list(c("CV","Area",2),c("Sum","SBR>5"),c("Sum","Height>0")),sType="St",addtoname="myST",addtoobject=TRUE)
head(mydata2$Annot)




}
\keyword{misc}
