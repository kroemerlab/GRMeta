\name{sortEntries}
\alias{sortEntries}

\title{Reorder metaboSet object entries}
\description{Reorder metaboSet object entries}
\usage{
sortEntries(obj, ...)
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

sortEntries(mydata,sType=TRUE,InjOrder=TRUE)$Meta
sortEntries(mydata,sType=FALSE,InjOrder=TRUE)$Meta

sortEntries(mydata,MetName=TRUE,IsSTD=FALSE)$Analyte
sortEntries(mydata,IsSTD=FALSE,MetName=TRUE)$Analyte


}

}
\keyword{misc}
