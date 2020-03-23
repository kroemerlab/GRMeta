loadSciexData<-function(ifile,...){
  
  #Sciex names of interest
  Samp.names.sc=c('Sample Name','Original Filename','Acquisition Date & Time')
  Mol.names.sc=c('Retention Time','Area','Height','Signal / Noise')
  
  #Agilent-Equivalent names
  Samp.names=c('Name','Data File','Acq. Date-Time')
  Mol.names=c('RT','Area','Height','S/N')
  
  #Read and reformat the table
  dm=read.table(ifile,sep='\t',header=T,check.names=F,stringsAsFactors=F)
  dl=lapply(unique(dm$`Sample Name`),function(x)dm[which(dm$'Sample Name'==x),]);names(dl)=unique(dm$`Sample Name`)
  Other.names.sc=setdiff(colnames(dm),c(Samp.names.sc,Mol.names.sc))
  
  #Turn format to Agilent-like
  am=do.call('rbind',lapply(dl,function(x){
    cbind(x[1,Samp.names.sc],do.call('cbind',plyr::alply(x[,c(Mol.names.sc,Other.names.sc)],1,function(y)y)))
  }))
  
  #Add column labels
  am=rbind(c(Samp.names,rep(c(Mol.names,Other.names.sc),length(dl[[1]]$`Component Name`))),am)
  
  #Add Molecule names
  header=c('Sample',rep('',length(Samp.names)-1),sapply(paste(dl[[1]]$`Component Name`,'Results'),function(x)c(x,rep('',length(c(Mol.names,Other.names.sc))-1))))
  am=rbind(header,am)
  
  #Export file to be read
  tp=tempfile()
  write.table(am,tp,row.names=F,col.names=F,quote = F,sep='\t')
  return(loadAgilentData(ifile=tp,...))
}
