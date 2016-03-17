\name{matchmzSet}
\alias{matchmzSet}
\title{Matching m/z to a database}
\description{
Matching m/z to a database
}
\usage{
matchmzSet(obj, Analyte = NULL, annotdb = NULL, ipdb = NULL, lIP = NULL, dppm = 20, mz2use = "MZ", mass2use = "Mass", 
annotinf2add = NULL,objinf2add = NULL, groupby =NULL, chunk = 200, collsep = ";;")

}
\arguments{
 \item{obj}{metaboSet or a data.frame}
  \item{Analyte}{Vector of analytes}
  \item{annotdb}{Dataframe containing the annotation infos}
  \item{ipdb}{Dataframe containing the IP definitions}
  \item{lIP}{List of IP, either from the Id, Name or Set columns}
  \item{dppm}{Delta ppm}
  \item{mz2use}{name of the column indicative of the m/z from the obj$Annot}
  \item{mass2use}{name of the column indicative of the mass in the annotdb}
  \item{annotinf2add}{column names in the annotdb that are added}
  \item{objinf2add}{column names in the  obj$Annot that are added}
  \item{groupby}{merging lines with the same values}
  \item{chunk}{split the obj in chunks: 100-500 is OK}
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

## explicit use of AnnotationDB
res=matchmzSet(mydata,lIP='[M+H]1+')

## all matches in IPDB$Set='Pos1'
res0=matchmzSet(mydata,lIP='Pos1')
print(table(res$IP))
print(table(res0$IP))

## add further infos
res1=matchmzSet(mydata,lIP='[M+H]1+',annotinf2add=c('GName','MF'))
## add further infos
res2=matchmzSet(mydata,lIP='[M+H]1+',annotinf2add=c('GName','MF'),objinf2add=c('RT','MZ'))
head(res1)
head(res2)


## collapse same 'MF'
res3=matchmzSet(mydata,lIP='[M+H]1+',annotinf2add=c('GName','MF'),groupby='MF')

head(res1)
head(res3)

## collapse on Analyte: not much sense if several IPs seek for!
res4=matchmzSet(mydata,lIP='[M+H]1+',annotinf2add=c('GName','MF'),groupby='Analyte')
table(table(res2$Analyte))
table(table(res4$Analyte))


}

}
\keyword{misc}
