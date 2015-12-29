\name{loadMavenData}
\alias{loadMavenData}
\title{Load MAVEN targetted data files}
\description{Load MAVEN targetted data files}
\usage{
loadMavenData(ifile, ofile = NULL,stdData=NULL, params = list())
}
\arguments{
  \item{ifile}{mzPeaks.mzroll file from Maven}
  \item{ofile}{Optional .Rdata file to store the object}
  \item{stdData}{Optional dataframe to calculate DRT and DPPM - requires GName, RT, IP fields as in \code{\link{stdDB}} }
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
testfile=list.files(system.file( package = "GRMeta"),pattern = "roll$",full.names = TRUE)

mydata=loadMavenData(ifile =testfile,params = list(AssayName="rp"))
print(mydata)

std=stdDB[stdDB$Method=="rppos",]
mydata=loadMavenData(ifile =testfile,stdData=std,params = list(AssayName="rp"))
print(mydata$Annot[1:10,1:7])

}

\keyword{misc}
