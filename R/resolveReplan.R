resolveReplan <- function(replan,B,N,M,k,import=NULL,export=NULL){
  ReName <- NULL
  k <- as.numeric(substring(grep("^cv[0-9]+$",replan,val=TRUE,ignore.case=TRUE),3))
  if (length(k)==0) k <- NULL
  if (is.null(k)){
    re.noinf <- length(grep("noinf",replan,val=FALSE,ignore.case=TRUE))>0
    re.boot <- length(grep("boot|plain",replan,val=FALSE,ignore.case=TRUE))>0
    re.outofbag <- length(grep("outofbag|out.of.bag|out-of-bag|b0|bootcv",replan,val=FALSE,ignore.case=TRUE))>0
    re.632 <- length(grep("632",replan,val=FALSE,ignore.case=TRUE))>0
    re.plus <- length(grep("plus|\\+",replan,val=FALSE,ignore.case=TRUE))>0
    ## replan <- match.arg(replan,c("none","plain","outofbag","boot632","boot.632","boot632plus","boot.632plus","noinf"))
    if (re.noinf==TRUE){replan <- "noinf"; ReName <- "no information error"}
    else if (re.outofbag==TRUE){replan <- "outofbag"; ReName <- "out-of-bag error"}
    else
      if (re.boot==TRUE){
        if (re.632==TRUE){
          if (re.plus==TRUE){replan <- "boot632plus"; ReName <- ".632+"}
          else{replan <- "boot632";ReName <- ".632"}
        }
        else{replan <- "plain"; ReName <- "bootstrap"}
      }
    if (is.null(ReName)) {replan <- "noPlan"; ReName <- "no plan"}
  }
  else{
    replan <- "crossval"
    ReName <- paste(k,"fold cross-validation",sep="-")
  }
  if (missing(M)) M <- N
  stopifnot(M>0 && M<=N) 
  subsampling <- M!=N
  ## if (missing(na.accept)) na.accept <- M/10
  if (replan=="noPlan"|| replan=="noinf") {
    B <- 0
  }
  else{
    if (missing(B)){
      if (length(k)>0) B <- 1 # repeat k-fold CrossVal ones
      else B <- 100  # either `plain' or `OutOfBag'
    }
    else if (B==0) stop("No. of resamples must be a positive integer.")
  }
  if (!is.null(import)){ # resampling data are provided by an external file
    if (is.list(import)){
      if (is.null(import$verbose)) import$verbose <- TRUE
      ResampleIndex <- do.call("pecImportIndex",import)
    }
    else{
      stopifnot(dim(import)==c(M,B))
      ResampleIndex <- import
    }
  }
  else{
    if (length(k)>0){
      ResampleIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
    }
    else{
      if (replan %in% c("boot632plus","outofbag","boot632","plain")){
        ResampleIndex <- matrix(sapply(1:B,function(b){sort(sample(1:N,size=M,replace=!subsampling))}),nrow=M,ncol=B)
        colnames(ResampleIndex) <- paste("Boot",1:B,sep=".")
      }
      else{
        ResampleIndex <- NULL
      }
    }
  }

  if (!is.null(ResampleIndex)){
    if (!is.null(export)){ # resampling data are exported to an external file
      stopifnot(is.list(export))
      if (is.null(export$verbose)) export$verbose <- TRUE
      do.call("pecExportIndex",c(list(x=ResampleIndex),export))
    }
  }
  out <- list(name=ReName,internal.name=replan,index=ResampleIndex,k=k,B=B,M=M,N=N)
  class(out) <- "rePlan"
  out
}

print.rePlan <- function(x){
  cat(paste("\nEstimation strategy:",x$name),"\n")
  if (x$internal.name=="crossval")
    cat("Repeat: ",x$B,"\n")
  else 
    cat("No. bootstrap samples: ",x$B,"\n")
  cat("Sample size: ",x$N,"\n")
}
