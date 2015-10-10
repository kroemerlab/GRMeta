\name{plot.metaboSet}
\alias{plot.metaboSet}
\title{Plot analyte data from a metaboSet object}
\description{Plot analytes matrices}
\usage{
\method{plot}{metaboSet}(obj,outfile=NULL,analyte=obj$Analyte[1],lgraphs=list(c("InjOrder","RT"),c("InjOrder","Area"),c("Height"),c("Height","Area")),
                         mfrow=c(2,2),deltaRT=0.05,linking="QC",orderBPlots="sType",...)
}
\arguments{
  \item{obj}{Experiment obj}
  \item{outfile}{PDF file: if not given results are graphed in the environment}
  \item{lgraphs}{List of size vectors for plotting x-y, if size 1 or the 2 elements are the same, a barlot is created}
  \item{mfrow}{Distribution of the plots on the figure}
  \item{deltaRT}{Censoring around the RT median}
  \item{linking}{Linking points from the specified sample type when injection order is on the x-axis}
  \item{orderBPlots}{Ordering barplots, default is "sType"}
  \item{\dots}{Arguments passed to \code{\link{par}} and \code{\link{pdf}}, excepted mfrow}
}
\value{Experiment object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{

testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)
mydata=loadAgilentData(ifile =testfiles[1])
print(mydata$Analyte)

analyte=mydata$Analyte[2]
\dontrun{
plot(obj,Analyte=mydata$Analyte[2],lgraphs=list(c("InjOrder","RT"),c("InjOrder","Area"),c("Height"),c("Height","Area")),
                         mfrow=c(2,2),deltaRT=0.05,linking="QC")
}
}
\keyword{misc}
