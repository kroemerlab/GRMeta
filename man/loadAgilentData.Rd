\name{loadAgilentData}
\alias{loadAgilentData}
\title{Load Agilent data text files}
\description{Load Agilent data text files}
\usage{
loadAgilentData(ifile, ofile = NULL, params = list())
}
\arguments{
  \item{ifile}{Tab separeted files}

  \item{ofile}{Optional .Rdata file to store the object}
  \item{params}{Parameters for loading/parsing and organising the data}
}
\value{Experiment data object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\seealso{
  \code{\link{paramsParsing}}
}
\examples{
testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)

mydata=loadAgilentData(ifile =testfiles[1])
print(mydata$Analyte)
print(mydata$Meta)

mydata=loadAgilentData(ifile =testfiles[1],params = list(AssayName="lc"))
print(mydata$Analyte)

mydata=loadAgilentData(ifile =testfiles[1],params = list(AssayName="lc",regTypes="^([BlStQC]+)_.*"))
print(mydata$Meta)
}

\keyword{misc}
