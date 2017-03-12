
######################################################################################
## Single point MZ calib over scan/mz

setGeneric(name=".GRcorrRawPPM",def=function(object,addppm) standardGeneric(".GRcorrRawPPM"))

setMethod(".GRcorrRawPPM","xcmsRaw",
          .GRcorrRawPPM<-function(object,addppm){
            ## single point calibration
            ## if vector, names must match the scan id
            
            #  xRaw2=xRaw
            scsten=cbind(object@scanindex+1,c(object@scanindex[-1],length(object@env$intensity)))
            rownames(scsten)=1:length(object@scanindex)
            
            ### Fix add ppm
            if(length(addppm)==1) addppm=unname(addppm)
            if(is.null(names(addppm))){
              addppm=rep(median(addppm),length(object@scanindex))
              names(addppm)=rownames(scsten)
            }
            
            newmz=object@env$mz
            for(ix in rownames(scsten)){
              oldy<-xcms:::getScan(object,as.numeric(ix))[,1]
              ppmy<- 1-addppm[ix]*10^-6
              #              ppmy<- 1-addppm*10^-6
              newy=round(oldy*ppmy,6)
              newmz[scsten[ix,1]:scsten[ix,2]]=newy
            }
            print(summary((1-newmz/object@env$mz)*10^6))
            
            ob<-new("xcmsRaw")
            ob@env <- new.env(parent = .GlobalEnv)
            ob@env$mz<-as.numeric(newmz)
            ob@env$intensity<-object@env$intensity
            ob@scanindex<-object@scanindex
            ob@scantime<-object@scantime
            
            ob@acquisitionNum<-1:length(ob@scanindex)
            ob@filepath<-object@filepath
            ob@mzrange<-range(ob@env$mz)
            ob@profmethod<-object@profmethod
            ob@tic<-object@tic
            ob@profparam<-list()
            ob<-xcms:::remakeTIC(ob)
            return(ob)
          }
          
)
######################################################################################

######################################################################################
## RT correction

setGeneric(name=".GRcorrRawRT",def=function(object,newrt) standardGeneric(".GRcorrRawRT"))

setMethod(".GRcorrRawRT","xcmsRaw",
          .GRcorrRawRT<-function(object,newrt){
            ## replace object@scantime with newrt + remove anything that is <=0 or is.na
            
            scsten=cbind(object@scanindex+1,c(object@scanindex[-1],length(object@env$intensity)))
            rownames(scsten)=1:length(object@scanindex)
            
            l2use=which(newrt>0 & !is.na(newrt))
            val2use=range(scsten[l2use,])
            scanindex=scsten[l2use,1]
            scanindex=scanindex-min(scanindex)
            
            ob<-new("xcmsRaw")
            ob@env <- new.env(parent = .GlobalEnv)
            ob@env$mz<-as.numeric(object@env$mz)[val2use[1]:val2use[2]]
            ob@env$intensity<-object@env$intensity[val2use[1]:val2use[2]]
            ob@scanindex<-as.integer(scanindex)
            ob@scantime<-newrt[l2use]
            
            ob@acquisitionNum<-1:length(scanindex)
            ob@filepath<-object@filepath
            ob@mzrange<-range(ob@env$mz)
            ob@profmethod<-object@profmethod
            ob@tic<-object@tic[l2use]
            ob@profparam<-list()
            ob<-xcms:::remakeTIC(ob)
            return(ob)
          }
          
)
