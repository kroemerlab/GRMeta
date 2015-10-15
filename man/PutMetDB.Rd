% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{PutMetDB}
\docType{data}
\alias{PutMetDB}

\title{Database of molecular formula}
\usage{
data(PutMetDB)
}

\description{
 Annotation database
}

\value{
Data frame:
 \item{GName}{Given name}
 \item{MF}{Molecular formula}
 \item{Mass}{Exact monoisotopic mass of neutral molecule}
 \item{HMDB}{HMDB ids}
 \item{KEGG}{Kegg ids}
 \item{CHEBI}{ChEbi ids}
}


\author{
David Enot \email{david.enot@gustaveroussy.fr}
}


\examples{
str(PutMetDB)
}    
\keyword{datasets}
