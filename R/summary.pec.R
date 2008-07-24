"summary.pec" <-  function(object,times,what,models,digits=3,print=TRUE,...){
  
  if (missing(models)) models <- names(object$models)
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),val=TRUE)
  }
  
  ##   NF <- 1:length(object)
  ##   if (!is.null(object$plainBoot)){
  ##     plainOut <- lapply(1:NF,function(f){
  ##       plainBoot.f <- object$plainBoot[[f]]
  ##       seBoot.f <- apply(plainBoot.f,1,function(b) sqrt(var(b,na.rm=TRUE)))
  ##       q.lower <- apply(plainBoot.f,1,function(b) quantile(b,0.025,na.rm=TRUE))
  ##       q.upper <- apply(plainBoot.f,1,function(b) quantile(b,0.975,na.rm=TRUE))
  ##       w.lower <- apparent.error - qnorm(.975) * seBoot.f
  ##       w.upper <- apparent.error + qnorm(.975) * seBoot.f
  ##     })
  ##   }
  
  cat("\nPrediction error curves\n\n")
  ##   cat(paste("\nEstimation method:",names(object$method),"\n"))
  print(object$method)
  otime <- object$time
  if (missing(times) && (length(times <- otime) > 20)){
    warning("Missing times argument: prediction error curves evaluated at the quantiles of fitted times\n")
    times <- quantile(otime)
  }
  tindex <- sindex(jump.times=object$time,eval.times=times)
  out <- lapply(what,function(w){
    cat("\n",w,"\n")
    tmp <- rbind(0, do.call("cbind",object[[w]][models]))[tindex+1,,drop=FALSE]
    tmp <- cbind(time=times,n.risk=c(object$n.risk[1],object$n.risk)[tindex+1],tmp)
    rownames(tmp) <- 1:NROW(tmp)
    if (print==TRUE) prmatrix(round(tmp,digits=digits),...)
  })
  cat("\n")
  invisible(out)
}
