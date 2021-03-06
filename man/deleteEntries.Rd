\name{deleteEntries}
\alias{deleteEntries}

\title{Remove entries from a metaboSet object}
\description{Remove entries from a metaboSet object}
\usage{
deleteEntries(obj, ...,verbose=TRUE)
}
\arguments{
  \item{obj}{metaboSet}
  \item{\dots}{see example}
  \item{verbose}{Show deleted entries?}
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
deleteEntries(mydata,MetName="CoA")$Annot

str(mydata$Data$Area)
newdata=deleteEntries(mydata,Sid=c("Bl_prerun-eau_1", "Bl_prerun-eau_1b", "Bl_prerun-matrice_4" ),MetName="Cytidine")
str(newdata$Data$Area)
}

}
\keyword{misc}
