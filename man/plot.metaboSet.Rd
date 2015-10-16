\name{plot.metaboSet}
\alias{plot.metaboSet}
\title{Plot analyte data from a metaboSet object}
\description{Plot analytes matrices}
\usage{
\method{plot}{metaboSet}(obj,outfile=NULL,analyte=obj$Analyte[1],
 lgraphs=list(c("RT~InjOrder"),c("Area~InjOrder",log="y"),c("Area~1",log="y"),c("Height~Area",log="xy")),
                         mfrow=c(2,2),deltaRT=0.05,linking="QC",orderBPlots=c("sType","Injorder"),cexBP=0.5,...)
}
\arguments{
  \item{obj}{Experiment obj}
  \item{outfile}{PDF file: if not given results are graphed in the environment}
  \item{lgraphs}{List of size vectors for plotting in form of y~x. Optionally 2 second argments is for controlling the log scale}
  \item{mfrow}{Distribution of the plots on the figure}
  \item{deltaRT}{Censoring around the RT median}
  \item{linking}{Linking points from the specified sample type when injection order is on the x-axis}
  \item{orderBPlots}{Ordering barplots, default is sType foollowed by injection order}
  \item{cexBP}{cex.axis value for labelling x-axis of the barplots}
  \item{\dots}{Arguments passed to \code{\link{par}} and \code{\link{pdf}}, excepted mfrow}
}
\value{Experiment object}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}

\examples{

testfiles=list.files(system.file( package = "GRMeta"),pattern = "txt$",full.names = TRUE)
mydata=loadAgilentData(ifile =testfiles[1])
mydata$Analyte

analyte=mydata$Analyte[2]
\dontrun{
plot(mydata,Analyte=mydata$Analyte[2],
    lgraphs=list(c("RT~InjOrder"),c("Area~InjOrder",log="y"),c("Height"),c("Height~Area",log="xy")),
                         mfrow=c(2,2),deltaRT=0.05,linking="QC")
                         
}
\dontrun{
#### sType: character
#### sTypeFac: factor
mydata$Meta$sTypeFac=factor(mydata$Meta$sType)
plot(mydata,Analyte=mydata$Analyte[2],
    lgraphs=list(c("Area",log="y"),c("Area~Sid",log="y"),
      c("Area~sType",log="y"),c("Area~sTypeFac",log="y")),
      mfrow=c(2,2),orderBPlots=c("sType","InjOrder"))
      
#### Only 1 boxplot -> mfrow=c(1,1)
plot(mydata,Analyte=mydata$Analyte[2],
    lgraphs=list(c("Area~sTypeFac",log="y")),
      mfrow=c(1,1))
}
#### Dump them all in a file
\dontrun{
plot(mydata,outfile="DumpAll.pdf",Analyte=mydata$Analyte,cex.main=1,height=8,cexBP = 0.4)
}


}
\keyword{misc}
