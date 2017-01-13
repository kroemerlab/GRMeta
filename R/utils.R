### check if x falls in the range or adjust to fit 
.GRchkrange<-function(x,rangeSc=c(1,Inf)){
  if(is.vector(x)){
    x=range(x,na.rm=T)
    return(c(max(floor(x[1]),rangeSc[1]),min(ceiling(x[2]),rangeSc[2])))
  }
  x=t(apply(x,1,range))
  x[,1]= floor(x[,1]) %>% {.[.<rangeSc[1]] = rangeSc[1]; .} 
  x[,2]= ceiling(x[,2]) %>% {.[.>rangeSc[2]] = rangeSc[2]; .} 
  x
}

# .GRchkrange2<-function(x,rangeSc=c(1,1,Inf,Inf)){
#   if(is.vector(x)){
#     x=c(floor(min(x,na.rm=T)),ceiling(max(x,na.rm=T)))
#     if(!(x[2]>rangeSc[2] | x[1]<rangeSc[3])) return(c(NA,NA))
#     ifelse(x[1])
#        
#     return(c(max(max(x[1],rangeSc[1]),rangeSc[2]),min(x[2]),rangeSc[2])))
#   }
#   x=t(apply(x,1,range))
#   x[,1]= floor(x[,1]) %>% {.[.<rangeSc[1]] = rangeSc[1]; .} 
#   x[,2]= ceiling(x[,2]) %>% {.[.>rangeSc[2]] = rangeSc[2]; .} 
#   x
# }

##### 
printElapse<-function(strtime,bef="",add="\n"){
  if(typeof(strtime)=="integer"){
    top=as.integer(as.POSIXct( Sys.time() ))-strtime
    units="secs"
  } else{
    top=Sys.time()-strtime
    units=attr(top,"units")
  }
  sprintf("%s%.2f %s%s",bef,top,units,add)
}

