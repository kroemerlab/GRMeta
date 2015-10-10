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
mydata1=loadAgilentData(ifile =testfiles[2])
mydata1
mydata2=loadAgilentData(ifile =testfiles[3])
mydata2

## merge 1 and 2
mydatacmb=mergeSet(mydata1,mydata2)
mydatacmb

## order matters
mydatacmb=mergeSet(mydata2,mydata1)
mydatacmb

######## all assays
mydata3=loadAgilentData(ifile =testfiles[1],params = list(AssayName="assay1"))
mydatacmb=mergeSet(mydata1,mydata2)
print(mydatacmb$Analyte)

}
\keyword{misc}
