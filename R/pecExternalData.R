pecExportIndex <- function(x,path,pattern,Fun,args,verbose=FALSE){
  if (!file.exists(path)) stop("Path in export list does not exist!")
  if (!(substring(path,nchar(path))=="/")) path <- paste(path,"/",sep="")
  if (missing(Fun) || is.null(Fun)) Fun <- "write.table"
  if (missing(args) || is.null(args)) args <- list(sep=";",row.names=FALSE,col.names=TRUE)
  if (missing(pattern) || is.null(pattern)) pattern <- "ResampleIndex"
  fname <- paste(path,pattern,sep="")
  do.call(Fun,c(list(x,file=fname),args))
  if (verbose==TRUE)
    message(paste("ResampleIndex saved to",
                  fname,
                  "\nusing",
                  Fun,
                  "\nwith extra arguments",
                  paste(names(args),args,collapse="=")))
}

pecImportIndex <- function(path,pattern,Fun,args,verbose=FALSE){
  if (!file.exists(path)) stop("Path in import list does not exist!")
  if (!(substring(path,nchar(path))=="/")) path <- paste(path,"/",sep="")
  if (missing(Fun) || is.null(Fun)) Fun <- "read.table"
  if (missing(args) || is.null(args)) args <- list(col.names=TRUE)
  if (missing(pattern) || is.null(pattern)) pattern <- "ResampleIndex"
  fname <- paste(path,pattern,sep="")
  if (!file.exists(fname)) stop(paste("File",fname,"does not exist"))
  x <- do.call(Fun,c(list(file=fname),args))
  if (verbose==TRUE)
    message(paste("ResampleIndex obtained from",
                  fname,
                  "\nusing",
                  Fun,
                  "\nwith extra arguments",
                  paste(names(args),args,collapse="="),"\nFound",NCOL(x),"bootstrap samples."))
  x
}
  
##   if (is.null(extern$data$pattern))
##     pattern <- c("TrainSet","TestSet")
##   else
##     pattern <- extern$data$pattern
##   if (is.null(extern$data$method)) extern$data$method <- "dir"
##   if (extern$data$method=="subdir"){
##     BootDir <- lapply(1:B,function(b){
##       dname <- paste(path,pecBootDataDir,sep=".")
##       stopifnot(dir.create(dname))
##       dname
##     })
##   }
##   else{
##     BootDir <- path
##   }
##   if (extern$data$create==TRUE){
##     if (is.null(extern$data$exportFun)) extern$data$exportFun <- "write.table"
##     if (is.null(extern$data$importFun)) extern$data$importFun <- "read.table"
##     if (is.null(extern$data$exportIndex) || extern$data$exportIndex){
##       bootFiles.extern <- FALSE
##       do.call(extern$data$exportFun,list(ResampleIndex,file=paste(BootDir,"ResampleIndex",sep="")))
##       message(paste("ResampleIndex saved to",paste(BootDir,"ResampleIndex",sep="")))
##     }
##     else{
##       bootFiles.extern <- TRUE
##       BootCVfiles <- lapply(1:B,function(b){
##         data$.internalId <- 1:N
##         vindex.b <- match(1:N,ResampleIndex[,b],nomatch=0)==0
##         BD <- if (extern$data$method=="subdir") BootDir[b] else BootDir
##         trainFile <- paste(BD,pattern[1],".",b,sep="")
##         testFile <- paste(BD,pattern[2],".",b,sep="")
##         val.b <- data[vindex.b,,drop=FALSE]
##         train.b <- data[ResampleIndex[,b],,drop=FALSE]
##         ##         do.call(extern$data$exportFun,list(train.b,file=trainFile))
##         ##         do.call(extern$data$exportFun,list(val.b,file=testFile))
##         list(trainFile=trainFile,testFile=testFile)
##       })
##     }
##   }
##   else{ ## read in existing bootstrap samples
##     stop("This functionality is not yet ready")
##     BootCVfiles <- lapply(1:B,function(b){
##       BD <- if (extern$data$method=="subdir") BootDir[b] else BootDir
##       trainFile <- paste(BootDir,pattern[1],".",b,sep="")
##       testFile <- paste(BootDir,pattern[2],".",b,sep="")
##       if (!(file.exists(trainFile))) stop(paste("Cant find file",trainFile))
##       if (!(file.exists(testFile))) stop(paste("Cant find file",testFile))
##       list(trainFile=trainFile,testFile=testFile)
##     })
##   }
## }
## }
