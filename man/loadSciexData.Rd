\name{loadSciexData}
\alias{loadSciexData}
\title{Load Sciex data text files}
\description{Load Sciex data text files}
\usage{
loadSciexData(ifile, ...)
}
\arguments{
  \item{ifile}{Tab separeted files}

  \item{...}{arguments to be passed to loadAgilentData() function}
}
\value{Experiment data object}
\author{
Allan Sauvat \email{allan.sauvat@gustaveroussy.fr}
}
\seealso{
  \code{\link{paramsParsing}}
}
\examples{
testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)

mydata=loadSciexData(ifile =testfiles[1])
print(mydata)

mydata=loadSciexData(ifile =testfiles[1],params = list(AssayName="lc"))
print(mydata)

}

\keyword{misc}
