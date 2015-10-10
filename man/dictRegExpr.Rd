% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{dictRegExpr}
\docType{data}
\alias{dictRegExpr}
\title{Regular expression dictionnary}
\description{Convenient dictionnary to ease match to the annotation database}
\usage{
data(dictRegExpr)
}

\value{
Data frame:
 \item{Before}{First argument of the \code{gsub} function}
 \item{After}{Second argument of the \code{gsub} function}
}

\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\seealso{
  \code{\link{cleanMetaboNames}}
}


\examples{
str(dictRegExpr)
}    
\keyword{datasets}
