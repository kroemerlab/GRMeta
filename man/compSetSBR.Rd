\name{compSetSBR}
\alias{compSetSBR}
\title{Add SBR matrix to the experiment data object}

\description{Compute signal to blank matrix and add to the metaboSet data object}

\usage{compSetSBR(obj, thresh = 1, SidBlank = NULL, what = "Height", newname = NULL, fct = median)}
\arguments{
  \item{obj}{Experiment object}
  \item{thresh}{Threshold used for NAs and censoring low intensities}
  \item{SidBlank}{Vector of Sids used for blank calculations}
  \item{what}{Matrix on which SBR calculation is performed}
  \item{newname}{Name of the SBR matrix to be added - if NULL returns the matrix rather than the object}
  \item{fct}{Function used to compute the blank level for each analyte}
}

\value{metaboSet object or list of SBR matrix and the normalisation factor}

\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{
(testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE))
mydata=loadAgilentData(ifile =testfiles[1])

mydata$Sid

newdata=compSetSBR(mydata,thresh=10,SidBlank=c( "Bl_prerun-eau_1","Bl_prerun-eau_1b", "Bl_prerun-matrice_4","Bl_postrun-matrice_4"),fct=median)
mydata

newdata=compSetSBR(mydata,thresh=10,SidBlank=NULL,newname="SBR")
newdata
head(newdata$Annot$SBR)

## return the SBR matrix  and not the metaboSet object
newdata=compSetSBR(mydata,thresh=10,SidBlank=NULL)
str(newdata)

}

\keyword{misc}
