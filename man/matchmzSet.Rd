\name{matchmzSet}
\alias{matchmzSet}
\title{Matching m/z to a database}
\description{
Matching m/z to a database
}
\usage{
matchmzSet(obj, Analyte = NULL, annotdb = NULL, ipdb = NULL, lIP = NULL, dppm = 20, mz2use = "MZ", mass2use = "Mass", infos2add = c("Gname", "MF"), groupby =NULL, chunk = 200, collsep = ";;")
}
\arguments{
 \item{obj}{metaboSet}
  \item{Analyte}{Vector of analytes}
  \item{annotdb}{Dataframe containing the annotation infos
}
  \item{ipdb}{Dataframe containing the IP definitions}
  \item{lIP}{List of IP, either from the Id, Name or Set columns}
  \item{dppm}{Delta ppm}
  \item{mz2use}{m/z from the obj$Annot}
  \item{mass2use}{mass column in the annotdb}
  \item{infos2add}{columns in the annotdb}
  \item{groupby}{merging lines with the same values}
  \item{chunk}{to be reduced for large obj}
  \item{collsep}{separator when merging lines}
}
\value{Data frame}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{
\dontrun{
testfile=list.files(system.file( package = "GRMeta"),pattern = "roll$",full.names = TRUE)

mydata=loadMavenData(ifile =testfile[1],params = list(AssayName="prof"))
print(mydata)

head(matchmzSet(mydata,annotdb=AnnotationDB,ipdb=IPDB,lIP='[M+H]1+'))

## add further infos
head(matchmzSet(mydata,annotdb=AnnotationDB,ipdb=IPDB,lIP='[M+H]1+',infos2add=c('GName','MF')))


## collapse same 'MF'
head(matchmzSet(mydata,annotdb=AnnotationDB,ipdb=IPDB,lIP='[M+H]1+',infos2add=c('GName','MF'),groupby='MF'))


}

}
\keyword{misc}
