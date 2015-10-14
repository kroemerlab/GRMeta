% wll-12-10-2006: First draft
% dle-drp-15-09-2007: Major changes

\name{stdDB}
\docType{data}
\alias{stdDB}

\title{Database of standards}
\usage{
data(stdDB)
}

\description{
 Annotation database
}

\value{
Data frame:
 \item{GName}{GName}
 \item{RT}{Retnetion time}
 \item{IP}{Mass of the ionisation product}
 \item{Method}{Assay id}
 \item{Type}{1: standard, 2: infered from MF, 3: to be tested}
 \item{MF}{Molecular formula}
 \item{MW}{Exact monoisotopic mass of neutral molecule}

}


\author{
David Enot \email{david.enot@gustaveroussy.fr}
}


\examples{
str(stdDB)
}    
\keyword{datasets}
