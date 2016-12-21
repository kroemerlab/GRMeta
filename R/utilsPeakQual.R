
.GRcodadw2<-function(x,lag=1) sum(diff(diff(x),lag=lag)^2)/sum(diff(x)^2)

.GRcodadw<-function(x) sum(diff(x)^2)/sum((x)^2)

### rewrite of ptw 1.9-11
.GRcoda<-function (x, window = 5, smoothing = c("median", "mean"),mcq=TRUE) 
{
  smoothing <- match.arg(smoothing)
  if (is.vector(x)) 
    x <- matrix(x, nrow = 1)
  x.smooth <- switch(smoothing, mean = t(apply(x, 1, function(xx) rowMeans(embed(xx, 
                                                                                 window)))), median = t(apply(x, 1, runmed, k = window, 
                                                                                                              endrule = "keep")))
  nc <- ncol(x)
  noff <- window%/%2
  if (smoothing == "median") 
    x.smooth <- x.smooth[, -c(1:noff, (nc - noff + 1):nc), 
                         drop = FALSE]
  x <- x[, -c(1:noff, (nc - noff + 1):nc), drop = FALSE]
  lambda <- sqrt(rowSums(x^2, na.rm = TRUE))
  A.lambda <- sweep(x, MARGIN = 1, STATS = lambda, FUN = "/")
  A.s <- t(scale(t(x.smooth), center = TRUE, scale = TRUE))
  if(mcq){
    mcq=rowSums(A.lambda * A.s)/sqrt(nc - window)
    if(is.na(mcq)) mcq=0
    return(mcq)
  } else invisible(A.s)
}
