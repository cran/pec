"print.pec" <- function(x,
                        times,
                        digits=3,
                        what=NULL,
                        ...){
  cat("\nPrediction error curves\n\n")
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
      cat("\nIPCW:",x$cens.model,"model")
    else cat("\nno censoring")}
  # }}}
  # {{{ discover what to print
  if (is.null(what))
    if (x$splitMethod$internal.name=="noPlan")
      what="AppErr"
    else
      what <- paste(x$splitMethod$internal.name,"Err",sep="")
  # }}}
  # {{{ echo estimation splitMethod
  if (!is.null(x$splitMethod)) print(x$splitMethod)
  if (missing(times)){
    ## times <- x$maxtime ## times <- quantile(x$time,.9)
    times <- x$minmaxtime
    ##     naPos <- sapply(x[[what]],function(pec){
    ##       length(pec)-sum(is.na(pec))-1
    ##     })
    ##     times <- min(x$time[naPos],times)
  }
  # }}}
  # {{{ cumulative prediction errors
  cat("\n",paste(rep("_",options()$width/2),collapse=""),"\n")
  tnames <- paste("time=",round(times,1),sep="")
  tnames[times<1] <- paste("time=",signif(times[times<1],2),sep="")
  cat("\nCumulative prediction error, aka Integrated Brier score  (IBS)\n aka Cumulative rank probability score\n\nRange of integration:",x$start,"and",tnames,":\n\n")
  out <- crps(object=x,times=times,start=x$start,what=what)
  if (is.matrix(out))
    print(out,digits=digits,quote=FALSE)
  else{
    print.listof(out,digits=digits,quote=FALSE)
  }
  # }}}
  # {{{ warn about failed computations
  failed <- sapply(x$failed,length)>0
  if (sum(failed)>0){
    cat("\nWarning: In some bootstrap samples the model failed:\n")
    print(sapply(x$failed[failed],table))
  }
  # }}}
  # {{{ echo test results
  if (!is.null(x$multiSplitTest)){
    cat("\n",paste(rep("_",options()$width/2),collapse=""),"\n")
    print(x$multiSplitTest)
  }
  # }}}
  invisible(x)
}
