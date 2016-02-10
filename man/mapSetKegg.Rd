\name{mapSetKegg}
\alias{mapSetKegg}
\title{Create KEGG map from object}
\description{Create KEGG map from object}
\usage{
mapSetKegg(obj,outfile=NULL, cols=NULL,sizes=4,mixcol="grey70",ftsize=16,conv2pdf=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{metaboSet object}
  \item{outfile}{File to save}
  \item{cols}{Colors used for each method in the object}
  \item{sizes}{Relative size for each node}
  \item{mixcol}{If null, mixing colors is done for entries in several methods}
  \item{ftsize}{Size of the label}
  \item{conv2pdf}{Convert SVG to PDF with Inkscape: command line}
}
\value{
SVG text
}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\examples{
\dontrun{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
lcdata=loadAgilentData(ifile =testfiles[1],params = list(AssayName="lc"))
gcdata1=loadAgilentData(ifile =testfiles[2],params = list(AssayName="gc"))

finaldata=mergeSet(lcdata,gcdata1)
mapSetKegg(finaldata,outfile="myKegg.svg")

}
}
\keyword{misc}
