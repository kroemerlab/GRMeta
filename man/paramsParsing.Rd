\name{paramsParsing}
\alias{paramsParsing}
\title{Default parameters for parsing Agilent text files}
\description{
Default parameters for parsing Agilent text files.
}
\usage{paramsParsing(AssayName = "myassay", FileCol = "Data File", TimeCol = "Acq. Date-Time", 
    ordering = TRUE, regTypes = "^([blBLQCcSTDstda]+)_.*", NameClean = c("_GCMRM", 
        "_MRM", "_DBAA"), checkNams = TRUE, Batch = NULL, nozeroscheck = NULL, 
    AnnotDB = AnnotationDB)}
\arguments{
 \item{AssayName}{Name of the method}
 \item{FileCol}{Name of the file column}
 \item{TimeCol}{Name of the file column}
 \item{ordering}{Reorder according to sample type and date of analysis}
 \item{regTypes}{Regular expression to extract sample type}
 \item{NameClean}{Strings being replaced from the filename in addition to .d and .D}
 \item{checkNams}{Verify analyte names (default call of \code{\link{cleanMetaboNames}})}
 \item{Batch}{Optional batch name for fileinfos}
 \item{nozeroscheck}{Optional vector of data matrices for which <=0 are not replace by NA. Default=NULL: all matrices will be checked}
 \item{AnnotDB}{Annotation DB: if null, load \code{\link{AnnotationDB}}}
}
\author{
David Enot \email{david.enot@gustaveroussy.fr}
}
\seealso{
  \code{\link{AnnotationDB}}
}

\examples{
names(paramsParsing())
}

