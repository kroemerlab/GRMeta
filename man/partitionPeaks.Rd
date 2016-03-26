\name{partitionPeaks}
\alias{partitionPeaks}
\title{Partition peak list}
\description{Partition peak list}
\usage{
partitionPeaks(cdat, err.ppm = 10, err.rt = 0.1, metaboSet = list(),refine=TRUE, fac = 4, whichrt = "rt", whichmz = "mz", whichsamp = "samp")
}
\arguments{
  \item{cdat}{data}
  \item{err.ppm}{ppm on mass}
  \item{err.rt}{rt error}
  \item{metaboSet}{list with Method, Meta, File}
  \item{refine}{refining step or loose split}
  \item{fac}{multiply the err.ppm/err.rt by a factor the loose step - speed up}
  \item{whichrt}{rt colname}
  \item{whichmz}{mz col name}
  \item{whichsamp}{sample name col}
}
\value{List or metaboSet}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{

}
\keyword{misc}
