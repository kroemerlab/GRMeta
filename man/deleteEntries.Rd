\name{deleteEntries}
\alias{deleteEntries}

\title{Remove entries from a metaboSet object}
\description{Remove entries from a metaboSet object}
\usage{
deleteEntries(obj, ...)
}
\arguments{
  \item{obj}{metaboSet}
  \item{\dots}{see example}
}
\value{New metaboSet}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{
\dontrun{
library(GRMeta)
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
mydata=loadAgilentData(ifile =testfiles[1])

deleteEntries(mydata,Sid=c( "St_ATP_2","St_ATP_2b","St_ATP_2c" ))$File
deleteEntries(mydata,PutativeName="CoA")$Annot

str(mydata$Data$Area)
newdata=deleteEntries(mydata,Sid=c("Bl_prerun-eau_1", "Bl_prerun-eau_1b", "Bl_prerun-matrice_4" ),PutativeName="Cytidine")
str(newdata$Data$Area)
}

}
\keyword{misc}
