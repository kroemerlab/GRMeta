% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{AnnotFluxGC}
\docType{data}
\alias{AnnotFluxGC}

\title{Annotation data for GC fluxomics}
\usage{
data(AnnotFluxGC)
}

\description{
 Annotation database
}

\value{
Data frame:
 \item{GName}{Given name}
 \item{MF}{Molecular formula of the precursor molecule}
 \item{Mass}{Exact monoisotopic mass of the precursor molecule}
 \item{Deriv}{Derivation procedure}
 \item{FragType}{Type of precursor molecule from the original derivatised molecule}
 \item{Q1th}{Expected Q1 from the method}
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
