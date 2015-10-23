\name{exportSetXLSX}
\alias{exportSetXLSX}
\title{Export object to Xlsx}
\description{Export object to Xlsx}
\usage{
exportSetXLSX(obj,outfile,ldata=names(obj$Data),sortSample=c("sType","InjOrder"),sortAnalyte=c("Method","RT"),transpose=FALSE,nround=3,characterNA=c("-","NA"))
}
\arguments{
  \item{obj}{metaboSet}
  \item{outfile}{Output file}
  \item{ldata}{Matrices from obj$Data}
  \item{sortSample}{Sorting sample based on obj$Meta names}
  \item{sortAnalyte}{Sorting analytes based on obj$Annot names}
  \item{transpose}{Transposing obj$Data matrices}
  \item{nround}{Rounding numerical data}
  \item{characterNA}{Replacing NAs by characterNA[1] for analyte/sample infos sheets and characterNA[2] for obj$Data matrices}
}
\value{Ordering of the samples and analytes}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{

#### best to increase the heap size before calling library(GRMeta)
options(java.parameters = "-Xmx8000m")

testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)

mydata=loadAgilentData(ifile =testfiles[1])
myratios=compSetRatios(mydata,l2add="Area",conv2set=TRUE)

\dontrun{
exportSetXLSX(mydata,outfile="test.xlsx",ldata=c("Area","RT"),characterNA=c("Missing","NAssss"))
exportSetXLSX(myratios,outfile="testratios.xlsx",ldata=c("Area","RT"),characterNA=c("Missing","NAssss"))
}

}

\keyword{misc}
