% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{dictSyno}
\docType{data}
\alias{dictSyno}
\title{Synonyms dictionnary}
\description{Convenient dictionnary to ease match to the annotation database}
\usage{
data(AnnotationDB)
}

\value{
Data frame:
 \item{Before}{Metabolite to be replaced}
 \item{After}{New metabolite name}
}

\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\seealso{
  \code{\link{cleanMetaboNames}}
}


\examples{
str(dictSyno)
}    
\keyword{datasets}
