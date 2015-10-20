\name{loadAgilentDataFlux}
\alias{loadAgilentDataFlux}
\title{Load Agilent data text files}
\description{Load Agilent data text files}
\usage{
loadAgilentDataFlux(ifile, ofile = NULL, params = list())
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

mydata=loadAgilentDataFlux(ifile =testfiles[4])
print(mydata)

}

\keyword{misc}
