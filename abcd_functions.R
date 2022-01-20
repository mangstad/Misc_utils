read.abcd = function(file,sep="\t",skip=1,cols=NULL,descriptions=FALSE,descriptions.only=FALSE) {
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
    if (descriptions.only) {
      temp = data[1,]
      temp[1,] = descrip
      temp
    } else {
      list(data=data,descrip=descrip)
    }
  } else {
    data
  }
}

multi.merge = function(...,by=NULL) {
  Reduce(function(x,y) merge(x,y, all=TRUE,by=by), list(...))
}

##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}
