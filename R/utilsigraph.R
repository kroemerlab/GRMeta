
## need igraph installed

## from gRbase v 1.8-1
.GRcombnPrim<-function (x, m, transpose = TRUE) 
{
  .colmat2list<-function (XX_) .Call("gRbase_colmat2list", PACKAGE = "GRMeta", XX_)
  if (length(x) == 1 && is.numeric(x)) 
    x <- seq(x)
  if (length(x) < m) 
    stop("Error in combnPrim: n < m\n")
  NCAND <- length(x)
  NSEL <- as.integer(m)
  NSET <- as.integer(choose(NCAND, NSEL))
  ANS <- rep.int(0L, NSET * NSEL)
  res <- .C("combnC", NSEL, NCAND, NSET, ANS, PACKAGE = "GRMeta")[[4]]
  #if (simplify) {
  res<-  matrix(x[res], nrow = NSEL, ncol = NSET)
  if(transpose) res<-t(res)
  res
  #}
  # else {
  #   res <- matrix(x[res], nrow = NSEL, ncol = NSET)
  #   res <- colmat2list(res)
  #   names(res) <- NULL
  #   res
  # }
}

.GRmergellx<-function(ll,nmin=2){
  
  ledges=do.call("rbind",lapply(ll,.GRcombnPrim,2))
  clu=igraph:::clusters(igraph:::graph_from_edgelist(ledges,directed=F))
  clu=tapply(1:length(clu$membership),clu$membership,c)
  clu[sapply(clu,length)>=nmin]
}

.GRisover<-function(vst,ven,retOne=FALSE,thr=0){
  isover=which(outer(ven,vst,"-")>thr & outer(vst,ven,"-")<thr,arr=T)
  if(all(isover[,1]==isover[,2])) clu=as.list(1:length(vst)) else{
    clu=igraph:::clusters(igraph:::graph_from_edgelist(isover,directed=F))
    clu=tapply(1:length(clu$membership),clu$membership,c)
  }
  clu[sapply(clu,length)>(!retOne)]
}

