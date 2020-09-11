read.abcd = function(file,sep="\t",skip=1,cols=NULL,descriptions=FALSE) {
  headers = names(read.table(file,sep=sep,header=T,stringsAsFactors=F)[-1,])
  if (descriptions) {
    descrip = names(read.table(file,sep=sep,header=T,stringsAsFactors=F,skip=1))
  }
  data = read.table(file,sep=sep,header=T,stringsAsFactors=F,skip=skip)
  names(data) = headers
  if (!is.null(cols)) {
    data = subset(data,select=cols)
  }
  if (descriptions) {
    list(data=data,descrip=descrip)
  } else {
    data
  }
}

multi.merge = function(...,by=NULL) {
  Reduce(function(x,y) merge(x,y, all=TRUE,by=by), list(...))
}
