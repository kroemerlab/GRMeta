\name{partitionPeaks}
\alias{partitionPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
partitionPeaks(cdat, err.ppm = 10, err.rt = 0.1, metaboSet = TRUE, fac = 4, whichrt = "rt", whichmz = "mz", whichsamp = "samp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdat}{
%%     ~~Describe \code{cdat} here~~
}
  \item{err.ppm}{
%%     ~~Describe \code{err.ppm} here~~
}
  \item{err.rt}{
%%     ~~Describe \code{err.rt} here~~
}
  \item{metaboSet}{
%%     ~~Describe \code{metaboSet} here~~
}
  \item{fac}{
%%     ~~Describe \code{fac} here~~
}
  \item{whichrt}{
%%     ~~Describe \code{whichrt} here~~
}
  \item{whichmz}{
%%     ~~Describe \code{whichmz} here~~
}
  \item{whichsamp}{
%%     ~~Describe \code{whichsamp} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cdat, err.ppm = 10, err.rt = 0.1, metaboSet = TRUE, 
    fac = 4, whichrt = "rt", whichmz = "mz", whichsamp = "samp") 
{
    system.time(lsp <- .GRgenpkl(cdat[, whichmz], cdat[, whichrt], 
        err.ppm = err.ppm * fac, err.rt = err.rt * fac * 2))
    cat("Loose split: num. peaks:", nrow(cdat), " -> num. clusters: ", 
        length(lsp), "\n", sep = "")
    system.time(lsp1 <- lapply(lsp, function(x) .GRsplitBoth(cdat[x, 
        whichmz], cdat[x, whichrt], x, bwrt = err.rt, bwppm = err.ppm)))
    lsp1 = unlist(lsp1, FALSE)
    lsp1 = lapply(lsp1, function(x) x[order(cdat[x, whichsamp], 
        cdat[x, whichsamp])])
    mrt = round(sapply(lsp1, function(x) median(cdat[x, whichrt])), 
        1)
    mrz = sapply(lsp1, function(x) median(cdat[x, whichmz]))
    ana = sprintf("\%.4f@\%.2f", mrz, mrt)
    if (any(table(ana) > 1)) {
        l = which(ana \%in\% names(which(table(ana) > 1)))
        ana[l] = sprintf("\%.5f@\%.3f", mrz, mrt)[l]
    }
    names(lsp1) = ana
    lso = order(mrt, mrz)
    lsp1 = lsp1[lso]
    cat("Density based split: num. peaks=", nrow(cdat), " -> num. clusters:", 
        length(lsp1), "\n", sep = "")
    if (!is.list(metaboSet) | is.null(metaboSet)) 
        return(lsp1)
    mrt = sapply(lsp1, function(x) median(cdat[x, whichrt]))
    mrz = sapply(lsp1, function(x) median(cdat[x, whichmz]))
    ana = sprintf("\%.4f@\%.2f-\%s", mrz, mrt, metaboSet$Method)
    if (any(table(ana) > 1)) {
        l = which(ana \%in\% names(which(table(ana) > 1)))
        ana[l] = sprintf("\%.5f@\%.3f-\%s", mrz, mrt, metaboSet$Method)[l]
    }
    names(lsp1) = ana
    Annot = data.frame(Analyte = ana, MetName = NA, IsSTD = FALSE, 
        LevelAnnot = 4, Method = metaboSet$Method, MZ = mrz, 
        RT = mrt, stringsAsFactors = FALSE)
    rownames(Annot) = ana
    lsids = unique(cdat[, whichsamp])
    Meta = File = data.frame(Sid = lsids, sType = "Sa", stringsAsFactors = F)
    Meta$InjOrder = 1:nrow(Meta)
    File$Date = chron(1:nrow(Meta))
    rownames(Meta) = rownames(File) = lsids
    if (!is.null(metaboSet$Meta)) 
        if (all(lsids \%in\% metaboSet$Meta$Sid)) 
            Meta = metaboSet$Meta[lsids, ]
    if (!is.null(metaboSet$File)) 
        if (all(lsids \%in\% metaboSet$File$Sid)) 
            File = metaboSet$File[lsids, ]
    mrt = sapply(lsp1, function(x) tapply(cdat[x, whichrt], cdat[x, 
        whichsamp], median, na.rm = T)[lsids])
    mmz = sapply(lsp1, function(x) tapply(cdat[x, whichmz], cdat[x, 
        whichsamp], median, na.rm = T)[lsids])
    dimnames(mrt) = dimnames(mmz) = list(lsids, names(lsp1))
    Data = list(RT = mrt, MZ = mmz)
    lvars = names(cdat)[sapply(cdat, class) \%in\% c("numeric", 
        "integer")]
    lvars = lvars[!lvars \%in\% c(whichrt, whichmz)]
    for (i in lvars) {
        m = sapply(lsp1, function(x) tapply(cdat[x, i], cdat[x, 
            whichsamp], sum, na.rm = T)[lsids])
        dimnames(m) = list(lsids, names(lsp1))
        Data[[i]] = m
    }
    allmat = list(Method = metaboSet$Method, Sid = lsids, Analyte = Annot$Analyte, 
        Annot = Annot, Meta = Meta, File = File, Data = Data)
    class(allmat) = append(class(allmat), "metaboSet")
    invisible(allmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
