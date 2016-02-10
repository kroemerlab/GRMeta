######

paramsParsing<-function(AssayName="myassay",FileCol="Data File",TimeCol="Acq. Date-Time",ordering=TRUE,
                        regTypes="^([blBLQCcSTDstda]+)_.*",NameClean=c("_GCMRM","_MRM","_DBAA"),
                        checkNams=TRUE,Batch=NULL,nozeroscheck=NULL,AnnotDB=AnnotationDB){
#   print(AnnotDB)
#   if(is.null(AnnotDB)){
#     AnnotationDB<-NULL
#     data(AnnotationDB,package='GRMeta')
#     print(str(AnnotationDB))
#     AnnotDB<-AnnotationDB
#     rm('AnnotationDB')
#   }
#   print(str(AnnotDB))
  list(AssayName="myassay",FileCol="Data File",TimeCol="Acq. Date-Time",ordering=TRUE,
       regTypes="^([blBLQCcSTDstda]+)_.*",NameClean=c("_GCMRM","_MRM","_DBAA"),
       checkNams=TRUE,Batch=NULL,nozeroscheck=NULL,AnnotDB=AnnotDB)
}

