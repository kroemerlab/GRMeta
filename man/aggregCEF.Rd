\name{aggregCEF}
\alias{aggregCEF}
\title{
Agreggating CEF files
}
\description{
Agreggating CEF files
}
\usage{
aggregCEF(lfiles, minrt = -Inf, maxrt = Inf, maxmz = 1001.11, minmz = 0.79, ncl = 1, verbose = TRUE)
}
\arguments{
  \item{lfiles}{Vector of file location - use the names attributes as Sample id}
  \item{minrt}{min RT}
  \item{maxrt}{max RT}
  \item{maxmz}{max MZ}
  \item{minmz}{min MZ}
  \item{ncl}{num of clusters for parallel}
  \item{verbose}{print each step}
}
\value{Data frame}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{

}
\keyword{misc}
