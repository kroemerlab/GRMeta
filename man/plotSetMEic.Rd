\name{plotSetMEic}
\alias{plotSetMEic}
\title{Plot series of Eics}
\description{Plot series of Eics}
\usage{
plotSetMEic(obj,WhichRT="rtcor",groupCol=NULL,colorCol=NULL,addfile="./",endfile="-Eic",repDots="-",cexEL=0.6,...)
}
\arguments{
  \item{obj}{metaboSet object}
  \item{WhichRT}{retention time to consider}
  \item{groupCol}{Grouping factor to plot Eics}
  \item{colorCol}{Column in obj$Meta to specify colors to be used}
  \item{addfile}{What to add before the filename}
  \item{endfile}{What to add to the end of the file name}
  \item{repDots}{Replacing dots}
  \item{cexEL}{cex for the Eic Legend}
  \item{\dots}{Arguments 1) used to match analytes from obj$Eic$File and 2) passed to \code{\link{par}} and \code{\link{pdf}} including mfrow}
}
\note{Experimental: to be used with profiling data from DEn }
\value{Experiment object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{
#### Dump them all in a file:  by default all of them
\dontrun{
plot(mydata,outfile="DumpAll.pdf",cex.main=1,height=8,cexBP = 0.4)
}
}

\keyword{hplot}
