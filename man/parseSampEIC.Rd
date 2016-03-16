\name{parseSampEIC}
\alias{parseSampEIC}
\alias{parseOneSampEIC}
\title{
Parse MS raw datafile to extract EICs
}
\description{
Parse MS raw datafile to extract EICs
}
\usage{
parseSampEIC(matfile, tabeic, corrt = NULL, npad = 3, stepmz = 1/1000, verbose = TRUE, ncl = 1, mzdata = FALSE, chunk = 200, local = FALSE)
}
\arguments{
  \item{matfile}{Information data frame, columns: Sid, unique sample id, In:input file, Out=where to store}
  \item{tabeic}{tabeic  data frame, must contain Id, mzmin, mzmax, }
  \item{corrt}{Matrix, row=Sid, column=retention time to adjust}
  \item{npad}{num. of additional scan on each side}
  \item{stepmz}{round m/z}
  \item{verbose}{should blahblah}
  \item{ncl}{num. of clusters - needs library(snowfall)}
  \item{mzdata}{is the input file a mzdata - needs library(mzR)}
  \item{chunk}{control the }
  \item{local}{for debugging}
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
