\name{cleanMetaboNames}
\alias{cleanMetaboNames}
\title{Clean analyte names}
\description{Clean analyte names}
\usage{
cleanMetaboNames(metnam, RegExpr = dictRegExpr, Syno = dictSyno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{metnam}{Vector of characters}
  \item{RegExpr}{Data frame used to sequentially rename metabolite names based on regular expression}
  \item{Syno}{Data frame used to sequentially rename metabolite names by exact matching}
}
\value{
Vector of characters
}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{
### Default dictionaries
cleanMetaboNames(metnam = c("Beta.glucose","ACETYLCOA POS","ACETYLCOA PS"),RegExpr=dictRegExpr,Syno=dictSyno)
cleanMetaboNames(metnam = c("Beta.glucose","ACETYLCOA POS","ACETYLCOA PS"))
### Do not replace
cleanMetaboNames(metnam = c("Beta.glucose","ACETYLCOA POS","ACETYLCOA PS"),RegExpr=NULL,Syno=NULL)
}
\keyword{misc}
