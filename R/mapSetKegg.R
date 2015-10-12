mapSetKegg<-function(obj,outfile=NULL, cols=NULL,sizes=4,mixcol="grey70",ftsize=16,conv2pdf=FALSE){
  
  
  if(length(cols)!=length(obj$Method))  cols=brewer.pal(8,"Set2")[1:length(obj$Method)]
  if(length(sizes)==0) sizes=4
  if(length(sizes)!=length(obj$Method)) sizes=rep(sizes[1],length(obj$Method))
  if(is.null(names(cols))) names(cols)=(obj$Method)
  if(is.null(names(sizes))) names(sizes)=(obj$Method)
  
  svgfile=list.files(system.file(package = "GRMeta"),pattern = "SVGstuff.rda$",full.names = TRUE)
  load(svgfile)
  
  llannot=tapply(obj$Annot$KEGG,obj$Annot$Method,function(x) unique(unlist(strsplit(na.omit(c(x)),"[;/]"))))
  
  #unlist(llannot)[!unlist(llannot)%in%keggids[,1]]
  
  keggids=keggids[keggids[,1]%in%unlist(llannot),]
  keggids$cols=keggids$assay=keggids$size=""
  for(i in names(llannot)){
    keggids$cols[keggids[,1]%in%llannot[[i]]]=paste(keggids$cols[keggids[,1]%in%llannot[[i]]],cols[i])
    keggids$assay[keggids[,1]%in%llannot[[i]]]=paste(keggids$assay[keggids[,1]%in%llannot[[i]]],i)
    keggids$size[keggids[,1]%in%llannot[[i]]]=paste(keggids$size[keggids[,1]%in%llannot[[i]]],sizes[i])
  }
  keggids$cols=gsub("^ ","",keggids$cols)
  keggids$assay=gsub("^ ","",keggids$assay)
  keggids$size=gsub("^ ","",keggids$size)
  
  if(!is.null(mixcol)){
    lcols=rowMeans(col2rgb(mixcol))/255
    mixcol=rgb(lcols[1],lcols[2],lcols[3])
  }
  
  for(i in grep(" ",keggids$cols)){
    if(is.null(mixcol)){lcols=rowMeans(col2rgb(sort(strsplit(keggids$cols[i]," ")[[1]])))/255
    keggids$cols[i]=rgb(lcols[1],lcols[2],lcols[3])} else{keggids$cols[i]=mixcol}
    
    keggids$assay[i]=paste(sort(strsplit(keggids$assay[i]," ")[[1]]),collapse=" ")
    keggids$size[i]=mean(as.numeric(strsplit(keggids$size[i]," ")[[1]]))
  }
  keggids$size=as.numeric(keggids$size)
  
  
  csvg=keggsvg;ltext=list();lnmatch=""
  
  i=1
  
  for(i in 1:nrow(keggids)){
    for(ix in grep(keggids[i,1],csvg)){
      old=csvg[[ix+1]]
      
      ## chnage color
      ireg=regexpr("fill:#[a-z0-9]+",old,perl=T)
      new=gsub(regmatches(old,ireg),paste("fill:",tolower(keggids$cols[i]),sep=""),old)
      ## chnage size
      ireg=regexpr(" cy=\"[0-9]+",new,perl=T)
      cy=as.numeric(gsub(" cy=\"","",regmatches(new,ireg)))
      ireg=regexpr(" cx=\"[0-9]+",new,perl=T)
      cx=as.numeric(gsub(" cx=\"","",regmatches(new,ireg)))
      ireg=regexpr(" ry=\"[0-9]+",new,perl=T)
      newry=as.numeric(gsub(" ry=\"","",regmatches(new,ireg)))
      newr=paste(" ry=\"",newry*keggids$size[i],sep="")
      new=gsub(regmatches(new,ireg),newr,new)
      ireg=regexpr(" rx=\"[0-9]+",new,perl=T)
      newrx=as.numeric(gsub(" rx=\"","",regmatches(new,ireg)))
      newr=paste(" rx=\"",newrx*keggids$size[i],sep="")
      new=gsub(regmatches(new,ireg),newr,new)
      
      csvg[[ix+1]]=new
      
      what=keggids$knam[i]
      ltext=c(ltext,list(paste("<g><text> <tspan y=\"",cy+newry*keggids$size[i]+ftsize,"\" x=\"",cx-newrx*keggids$size[i],
                               "\" style=\"font-size:",ftsize,"px;font-style:normal;font-weight:normal;fill:",
                               tolower(keggids$col[i]),";font-family:Arial\">",what,"</tspan></text></g>",sep="")))
      
    }
  }
  
  csvg2=c(unlist(csvg[-length(csvg)]),unlist(ltext),csvg[[length(csvg)]])
  if(!is.null(outfile)){
    cat(csvg2,file=outfile,sep="\n")
    if(conv2pdf) system(paste("inkscape -A ",gsub(".[sS][vV][gG]$",".pdf$",outfile)," ",outfile))
    
  }
  invisible(csvg2)
  
}


