\name{exportDataXLSX}
\alias{exportDataXLSX}
\title{Export object to Xlsx}
\description{Export object to Xlsx}
\usage{
exportDataXLSX(obj,outfile,ldata=names(obj$Data),sortSid=c("sType","InjOrder"),sortAnalyte=c("Method","RT"),transpose=FALSE,nround=3,characterNA=c("-","NA"))
}
\arguments{
  \item{obj}{metaboSet}
  \item{outfile}{Output file}
  \item{ldata}{Matrices from obj$Data}
  \item{sortSid}{Sorting sample based on obj$Meta names}
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

testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)

mydata=loadAgilentData(ifile =testfiles[1])
\dontrun{
exportDataXLSX(mydata,outfile="test.xlsx",ldata=c("Area","RT"),characterNA=c("Missing","NAssss"))
}

}

\keyword{misc}
