\name{corrSetQC}
\alias{corrSetQC}
\title{Compute QC adjusted matrices}
\description{Compute QC adjusted matrices}
\usage{
corrSetQC(obj,what, Samp2Corr = obj$Sid, Var2Corr = obj$Analyte, lQC = obj$Sid[which(obj$Meta$sType == "QC")], nminQC = 3, propNNA = 0.5, lPcs = 1:2, outfile = NULL, doplot = TRUE, complete = "nothing", ipc = 1, imod = 4)
}
\arguments{
  \item{obj}{metaboSet object}
  \item{what}{Which matrix to use}
  \item{Samp2Corr}{List of Sids to be corrected}
  \item{Var2Corr}{List of Analyte used for correction}
  \item{lQC}{List of QC Sids}
  \item{nminQC}{Number of non-missing QC for each variables. Default: 3}
  \item{propNNA}{Proportion of non-missing allowed in the QC samples. Default: 0.5}
  \item{lPcs}{Vector of PCs to be investigated. Default: 1:2}
  \item{outfile}{Optional PDF file to dum p the plots}
  \item{doplot}{Do plots?}
  \item{complete}{If any Sid not in Samp2Corr or Analyte in Var2Corr excluded or missing. 
  nothing:  former data are used in place of the new data, NA: new data replaced by NA, remove: corresponding rows/columns are excluded}
  \item{ipc}{PC used for adjusting}
  \item{imod}{Model  used for adjusting: 1, linear 2, quadratic 3, order 3 4, spline}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Matrix if ipc and imod specified or the prcomp object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
## TO DO
}
\keyword{misc}
