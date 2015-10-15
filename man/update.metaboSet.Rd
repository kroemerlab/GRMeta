\name{update.metaboSet}
\alias{update.metaboSet}

\title{Update ids in a metaboSet object}
\description{Update sample or analyte ids in an experiment object}
\usage{
\method{update}{metaboSet}(obj,what="Sid",formerid=NULL,newid=NULL,exact=TRUE,swap=FALSE)
}
\arguments{
  \item{obj}{Metaboset object}
  \item{what}{Only Sid supported}
  \item{formerid}{Vector of current Sids}
  \item{newid}{Vector of new Sids}
  \item{exact}{Exact matching}
  \item{swap}{Should names in oldname be swaped with corresponding newname}
}
\details{
If necessary, more details than the description above ~~
}
\value{Experiment object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{

testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)
mydata=loadAgilentData(ifile =testfiles[1])
print(mydata$Sid)

### Exact=TRUE: simple renaming of the Sid
mydata2=update(mydata,what="Sid",formerid=c("Bl_prerun-eau_1b","Bl_prerun-matrice_4"),newid = c("Bl-pre","BlM-pre"),exact=TRUE,swap = FALSE)
print(mydata2$File)

### Exact=FALSE: replace patterns by the second one
mydata2=update(mydata,what="Sid",formerid=c("pre","eau"),newid = c("*****",""),exact=FALSE,swap = FALSE)
print(mydata2$Meta)

### Swap two entries does ONLY exchange ids but not the file/date etc...
print(mydata$File)
mydata2=update(mydata,what="Sid",formerid="Bl_prerun-eau_1b",newid = "Bl_prerun-matrice_4",swap = TRUE)
print(mydata2$File)

}
\keyword{misc}
