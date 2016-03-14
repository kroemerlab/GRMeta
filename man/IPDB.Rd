% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{IPDB}
\docType{data}
\alias{IPDB}

\title{Annotation database}
\usage{
data(IPDB)
}

\description{
 IP database
}

\value{
Data frame:
 \item{Id}{Short id}
 \item{Name}{Name}
 \item{Set}{Convenience set groups}
 \item{Charge}{Total charge}
 \item{xM}{Num. of molecules}
 \item{adMass}{Mass to add}
 \item{adAt}{Atoms to add}
 \item{rmAt}{Atoms to remove}
}


\author{
David Enot \email{david.enot@gustaveroussy.fr}
}


\examples{
str(IPDB)
}    
\keyword{datasets}
