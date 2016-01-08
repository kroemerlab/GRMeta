% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{AnnotFluxGC}
\docType{data}
\alias{AnnotFluxGC}
\alias{AnnotFluxGC2}

\title{Annotation data for GC fluxomics}
\usage{
data(AnnotFluxGC)
data(AnnotFluxGC2)
}

\description{
 Annotation database
}

\value{
Data frame:
 \item{GName}{Given name}
 \item{MF}{Molecular formula of derivatised molecule}
 \item{Mass}{Exact monoisotopic mass of derivatised molecule}
 \item{FragType}{Type of fragment on Q1}
 \item{Q1th}{Method Q1}
 \item{MFori}{Molecular formula of the original molecule}
 \item{Class}{Chemical class}
 \item{HMDB}{HMDB ids}
 \item{KEGG}{Kegg ids}
 \item{CHEBI}{ChEbi ids}
}


\author{
David Enot \email{david.enot@gustaveroussy.fr}
}


\examples{
str(AnnotFluxGC)
}    
\keyword{datasets}
