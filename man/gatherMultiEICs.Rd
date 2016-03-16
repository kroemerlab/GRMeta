\name{gatherMultiEICs}
\alias{gatherMultiEICs}
\alias{gatherOneEIC}
\title{Fetch EICs from multiple sample files into one file}
\description{
Fetch EICs from multiple sample files into one file
}
\usage{
gatherMultiEICs( matfile, tabeic, root = "./", outfile = NA, ncl = 1, nsplit = 20, serial = TRUE)
}
\arguments{
  \item{matfile}{Information data frame, columns: Sid, unique sample id, In:input file, Out=where to store}
  \item{tabeic}{tabeic  data frame, must contain Id, mzmin, mzmax, }
  \item{root}{Directory to write down}
  \item{outfile}{onot really used}
  \item{ncl}{num. of clusters - needs snowfall}
  \item{nsplit}{how many to do in the same time}
  \item{serial}{useless??}
}
\details{
blahbhla
}
\value{Data frame}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
print('none')
}
\keyword{misc}
