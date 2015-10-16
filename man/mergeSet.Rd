\name{mergeSet}
\alias{mergeSet}
\title{Merge data from the same method or  different methods}
\description{Merge data from the same method or  different methods}
\usage{
mergeSet(\dots)
}
\arguments{
  \item{\dots}{metaboSet objects}
}
\value{metaboSet where all experiments are collated}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
## Better to give an explicit name for the batch
mydata1=loadAgilentData(ifile =testfiles[2],params = list(AssayName="A1",Batch="B1"))
mydata1
mydata2=loadAgilentData(ifile =testfiles[3],params = list(AssayName="A1",Batch="B2"))
mydata2

## merge 1 and 2
mydatacmb=mergeSet(mydata1,mydata2)
str(mydatacmb$Meta)

## order matters
mydatacmb=mergeSet(mydata2,mydata1)
str(mydatacmb$Meta)

## Infos in Meta/Annot/File will carried forward
mydata1$Meta$Grp=factor(rbinom(nrow(mydata1$Meta),1,.5))
mydatacmb=mergeSet(mydata1,mydata2)
mydatacmb$Meta$Grp

## remerging -> duplicating same samples!!!
mydatacmb2=mergeSet(mydata2,mydatacmb)
(mydatacmb2$File)
mydatacmb2$Sid[order(mydatacmb2$Meta$InjOrder)]

######## two different assays
mydata3=loadAgilentData(ifile =testfiles[1],params = list(AssayName="A2"))
mydatacmb=mergeSet(mydata1,mydata3)
mydatacmb

######## alltogether
mydatacmb=mergeSet(mydata1,mydata2,mydata3)
mydatacmb

}
\keyword{misc}
