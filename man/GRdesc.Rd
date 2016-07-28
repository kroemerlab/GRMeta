% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{GRdesc}
\docType{data}
\alias{GRdesc}

\title{Annotation database descriptor}
\usage{
data(GRdesc)
}

\description{
 Annotation database
}

\value{
Data frame:
 \item{GName}{GName}
 \item{Class}{Chemical class}
 \item{MF}{Molecular formula}
 \item{Mass}{Exact monoisotopic mass}
 \item{HMDB}{HMDB}
 \item{KEGG}{KEGG}
 \item{PC}{PC}
}


\author{
David Enot \email{david.enot@gustaveroussy.fr}
}


\examples{
str(GRdesc)
}    
\keyword{datasets}
