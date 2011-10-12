print.Cindex <- function(x,
                         digits=3,
                         what=NULL,
                         times,
                         ...){
  cat("\nThe c-index for right censored event times\n\n")
  # {{{ echo models
  cat("Prediction models:\n\n")
  printModels <- sapply(x$models,function(m){
    if (class(m) %in% c("character","call"))
      m
    else
      if (class(try(m$call,silent=TRUE))=="try-error")
        "unknown formula"
      else
        m$call
  })
  print(printModels,quote=FALSE)
  # }}}
  # {{{ echo response
  print(x$response)
  # }}}
  # {{{ echo cens model
  if (!is.null(x$cens.model)){
    if (x$cens.model!="none")
      cat("\nCensoring model for IPCW:",x$cens.model,"model",ifelse(x$cens.model=="marginal","(Kaplan-Meier for censoring distribution)",""),"\n")
    else cat("\nno censoring")}
  if (!is.null(x$splitMethod)) print(x$splitMethod)
  # }}}
  # {{{ discover what to print
  if (missing(what) || is.null(what)){
    what <- grep(c("Cindex$"),names(x),val=TRUE)
  }
  # }}}
  # {{{ result table
  out <- lapply(what,function(r){
    out <- do.call("rbind",lapply(1:length(x$models),function(m){
      x[[r]][[m]]
    }))
    if (is.matrix(out)){
      rownames(out) <- names(x$models)
      coln <- paste("time=",round(x$time,1),sep="")
      coln[x$time<1] <- paste("time=",round(x$time[x$time<1],4),sep="")
      colnames(out) <- coln
    }
    out
  })
  names(out) <- what
  # {{{ if only one time point
  if (NCOL(out[[1]])==1){
    cat("\nEstimated C-index in % at",colnames(out[[1]]),"\n\n")
    outMat <- 100*do.call("cbind",out)
    colnames(outMat) <- what
    if (!is.null(x$Pairs))
      outMat <- cbind(outMat,Pairs=round(x$Pairs[[1]],1),Concordant=round(unlist(x$Concordant),1))
    print(outMat,digits)
  }
  # }}}
  # {{{ multiple time points
  else{
    cat("\nEstimated C-index in %\n\n")
    print(lapply(out,function(x)x*100),digits)
  }
  if(x$splitMethod$name=="BootCv")
    cat("\nAppCindex    : Apparent performance\nBootCvCindex : Bootstrap crossvalidated performance\n\n")
  # }}}
  invisible(out)
}
